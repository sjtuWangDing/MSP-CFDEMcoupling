/*---------------------------------------------------------------------------*\
  CFDEMcoupling - Open Source CFD-DEM coupling

  CFDEMcoupling is part of the CFDEMproject
  www.cfdem.com
                              Christoph Goniva, christoph.goniva@cfdem.com
                              Copyright 2009-2012 JKU Linz
                              Copyright 2012-     DCS Computing GmbH, Linz
------------------------------------------------------------------------------
License
  This file is part of CFDEMcoupling.

  CFDEMcoupling is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 3 of the License, or (at your
  option) any later version.

  CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with CFDEMcoupling; if not, write to the Free Software Foundation,
  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
  cfdemCloud class managing DEM data for CFD-DEM coupling

Class
  Foam::cfdemCloud
\*---------------------------------------------------------------------------*/

#include "cloud/cfdem_cloud.h"
#include "cloud/coupling_properties.h"
#include "cloud/particle_cloud.h"
#include "mpi.h"
#include "sub_model/averaging_model/averaging_model.h"
#include "sub_model/data_exchange_model/data_exchange_model.h"
#include "sub_model/force_model/force_model.h"
#include "sub_model/force_model/global_force.h"
#include "sub_model/liggghts_command_model/liggghts_command_model.h"
#include "sub_model/locate_model/locate_model.h"
#include "sub_model/mom_couple_model/implicit_couple.h"
#include "sub_model/mom_couple_model/mom_couple_model.h"
#include "sub_model/void_fraction_model/void_fraction_model.h"

namespace Foam {

cfdemDefineTypeName(cfdemCloud);

cfdemCloud::cfdemCloud(const fvMesh& mesh)
    : mesh_(mesh),
      couplingPropertiesDict_(IOobject("couplingProperties",  // coupling properties file name
                                       mesh.time().constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)),
      liggghtsCommandsDict_(IOobject("liggghtsCommands",  // liggghts commands file name
                                     mesh.time().constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)),
      cProps_(mesh, couplingPropertiesDict_, liggghtsCommandsDict_),
      parCloud_(0),
      writeTimePassed_(false),
      meshHasUpdated_(false),
      validCouplingStep_(false),
      globalForce_(globalForce::New(*this, couplingPropertiesDict_)),
      dataExchangeModel_(dataExchangeModel::New(*this, couplingPropertiesDict_)),
      voidFractionModel_(voidFractionModel::New(*this, couplingPropertiesDict_)),
      locateModel_(locateModel::New(*this, couplingPropertiesDict_)),
      averagingModel_(averagingModel::New(*this, couplingPropertiesDict_)),
#if defined(version24Dev)
      turbulence_(mesh.lookupObject<turbulenceModel>(cProps_.turbulenceModelType())),
#elif defined(version21) || defined(version16ext)
#ifdef compre
      turbulence_(mesh.lookupObject<compressible::turbulenceModel>(cProps_.turbulenceModelType())),
#else
      turbulence_(mesh.lookupObject<incompressible::turbulenceModel>(cProps_.turbulenceModelType())),
#endif
#elif defined(version15)
      turbulence_(mesh.lookupObject<incompressible::RASModel>(cProps_.turbulenceModelType())),
#endif
      turbulenceMultiphase_(
          IOobject("turbulenceMultiphase", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
#ifdef compre
          dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0)  // kg/m/s
#else
          dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0)  // m²/s
#endif
          ),
      ddtVoidFraction_(
          IOobject("ddtVoidFraction", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
          dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0), 0)) {
  // create liggghts command model
  for (const auto& name : liggghtsCommandModelList()) {
    // liggghtsCommandModel::New() 函数返回的是 std::unique_ptr(xvalue)
    liggghtsCommandModels_.emplace_back(liggghtsCommandModel::New(*this, liggghtsCommandsDict_, name));
  }
  // create force model
  for (const auto& name : forceModelList()) {
    forceModels_.emplace_back(forceModel::New(*this, couplingPropertiesDict_, name));
  }
  // create mom couple model
  for (const auto& name : momCoupleModelList()) {
    std::shared_ptr<momCoupleModel> sPtr(momCoupleModel::New(*this, couplingPropertiesDict_, name));
    momCoupleModels_.insert(std::make_pair(name, sPtr));
  }
  // check periodic
  if (checkPeriodicCells() != checkSimulationFullyPeriodic()) {
    FatalError << "checkSimulationFullyPeriodic(): " << (checkSimulationFullyPeriodic() ? true : false)
               << ", but from dictionary read checkPeriodicCells: " << (checkPeriodicCells() ? true : false)
               << abort(FatalError);
  }
  // 执行初始化函数
  init();
}

cfdemCloud::~cfdemCloud() {}

void cfdemCloud::init() {
  // 由于在构造 dataExchangeModel 的时候会执行 liggghts 脚本，所以此时可以获得最初的颗粒信息
  int number = dataExchangeM().getNumberOfParticlesFromDEM();
  setNumberOfParticles(number);
  // allocate memory of init data exchanged with liggghts
  dataExchangeM().realloc(parCloud_.initVelocities(), base::makeShape2(number, 3), parCloud_.initVelocitiesPtr(), 0.0);
  dataExchangeM().getData("v", "vector-atom", parCloud_.initVelocitiesPtr());
}

//! \brief 重新分配内存
void cfdemCloud::reallocate() {
  int number = numberOfParticles();
  // allocate memory of data exchanged with liggghts
  dataExchangeM().realloc(parCloud_.radii(), base::makeShape1(number), parCloud_.radiiPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.cds(), base::makeShape1(number), parCloud_.cdsPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.positions(), base::makeShape2(number, 3), parCloud_.positionsPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.velocities(), base::makeShape2(number, 3), parCloud_.velocitiesPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.DEMForces(), base::makeShape2(number, 3), parCloud_.DEMForcesPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.fluidVel(), base::makeShape2(number, 3), parCloud_.fluidVelPtr(), 0.0);
  // allocate memory of data not exchanged with liggghts
  parCloud_.particleOverMeshNumber() = std::move(base::CITensor1(base::makeShape1(number), 0));
  parCloud_.findCellIDs() = std::move(base::CITensor1(base::makeShape1(number), -1));
  parCloud_.dimensionRatios() = std::move(base::CDTensor1(base::makeShape1(number), -1.0));
  parCloud_.impForces() = std::move(base::CDTensor2(base::makeShape2(number, 3), 0.0));
  parCloud_.expForces() = std::move(base::CDTensor2(base::makeShape2(number, 3), 0.0));
  parCloud_.particleRootProcIDs() = std::move(base::CITensor1(base::makeShape1(number), -1));
}

void cfdemCloud::getDEMData() {
  dataExchangeM().getData("radius", "scalar-atom", parCloud_.radiiPtr());
  dataExchangeM().getData("x", "vector-atom", parCloud_.positionsPtr());
  dataExchangeM().getData("v", "vector-atom", parCloud_.velocitiesPtr());
}

void cfdemCloud::giveDEMData() const {
  dataExchangeM().giveData("dragforce", "vector-atom", DEMForcesPtr());
  bool implDEMdrag = false;
  for (const auto& ptr : forceModels_) {
    if (ptr->forceSubM()->treatDEMForceImplicit()) {
      implDEMdrag = true;
      break;
    }
  }
  if (implDEMdrag) {
    dataExchangeM().giveData("Ksl", "scalar-atom", cdsPtr());
    dataExchangeM().giveData("uf", "vector-atom", fluidVelPtr());
  }
}

void cfdemCloud::printParticleInfo() const {
  int nProcs = 0, id = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  base::MPI_Barrier();
  if (0 == id) {
    for (int index = 0; index < numberOfParticles(); ++index) {
      Pout << "  position[" << index << "]: " << positions()[index][0] << ", " << positions()[index][1] << ", "
           << positions()[index][2] << endl;
    }
    for (int index = 0; index < numberOfParticles(); ++index) {
      Pout << "  velocity[" << index << "]: " << velocities()[index][0] << ", " << velocities()[index][1] << ", "
           << velocities()[index][2] << endl;
    }
  }
  base::MPI_Barrier();
  for (int index = 0; index < numberOfParticles(); ++index) {
    Pout << "  dimensionRatio[" << index << "]: " << dimensionRatios()[index] << endl;
  }
  base::MPI_Barrier();
  for (int index = 0; index < numberOfParticles(); ++index) {
    Pout << "  DEMForce[" << index << "]: " << DEMForces()[index][0] << ", " << DEMForces()[index][1] << ", "
         << DEMForces()[index][2] << endl;
  }
  base::MPI_Barrier();
  for (int index = 0; index < numberOfParticles(); ++index) {
    Pout << "  impForce[" << index << "]: " << impForces()[index][0] << ", " << impForces()[index][1] << ", "
         << impForces()[index][2] << endl;
  }
  base::MPI_Barrier();
  voidFractionM().printVoidFractionInfo();
}

//! \brief reset field
void cfdemCloud::resetField() {
  // 重置局部平均颗粒速度
  averagingM().resetUs();
  Info << "Reset Us fields - done" << endl;

  // 重置颗粒速度影响因数场
  averagingM().resetUsWeightField();
  Info << "Reset Us weight fields - done" << endl;

  // 重置空隙率场
  voidFractionM().resetVoidFraction();
  Info << "Reset voidfraction fields - done" << endl;

  // 重置隐式力场
  globalF().resetImpParticleForce();
  Info << "Reset implicit force fields - done" << endl;

  // 重置显式力场
  globalF().resetExpParticleForcee();
  Info << "Reset Explicit force fields - done" << endl;

  // 重置单位体积动量交换场
  for (const auto& pair : momCoupleModels_) {
    (pair.second)->resetMomSourceField();
  }
  Info << "Reset Ksl fields - done" << endl;
}

/*!
 * \brief 更新函数
 * \note used for cfdemSolverPiso
 * \param U      <[in] 流体速度场
 * \param voidF  <[in, out] 小颗粒空隙率场
 * \param Us     <[in, out] 局部平均小颗粒速度场
 * \param Ksl    <[in, out] 动量交换场
 */
void cfdemCloud::evolve(volVectorField& U, volScalarField& voidF, volVectorField& Us, volScalarField& Ksl) {
  Info << "/ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /" << endl;
  // 检查当前流体时间步是否同时也是耦合时间步
  validCouplingStep_ = dataExchangeM().checkValidCouplingStep();
  if (validCouplingStep_) {
    // 创建用于记录 coupling time step counter
    auto pCounter = std::make_shared<dataExchangeModel::CouplingStepCounter>(dataExchangeM());
    // couple(): run liggghts command and get number of particle
    setNumberOfParticles(dataExchangeM().couple());
    // reset field
    resetField();
    // realloc memory
    reallocate();
    // 获取 DEM data
    getDEMData();
    // 获取到在当前 processor 上颗粒覆盖的某一个网格编号，如果获取到的网格编号为 -1，则表示颗粒不覆盖当前 processor
    locateM().findCell(parCloud_.findCellIDs());
    // 计算颗粒空隙率
    voidFractionM().setVoidFraction();
    voidF = voidFractionM().voidFractionInterp();
    // 计算局部平局颗粒速度场
    averagingM().setVectorFieldAverage(averagingM().UsNext(), averagingM().UsWeightField(), velocities(),
                                       particleWeights());
    Us = averagingM().UsInterp();
    // global force init
    globalF().initBeforeSetForce();
    // 计算流体对颗粒的作用力
    for (const auto& ptr : forceModels_) {
      ptr->setForce();
    }
    // global force end
    globalF().endAfterSetForce();
    // 计算局部累加的流体作用力场
    averagingM().setVectorFieldSum(globalF().impParticleForce(), impForces(), particleWeights());
    // write DEM data
    giveDEMData();
    // get shared ptr of implicitCouple model and update Ksl field
    std::shared_ptr<momCoupleModel> sPtr = momCoupleModels_[implicitCouple::cTypeName()];
    Ksl = sPtr->impMomSource();
    Ksl.correctBoundaryConditions();
    printParticleInfo();
  }
  Info << __func__ << " - done\n" << endl;
}

/*!
 * \brief check if simulation is fully periodic
 * \return true if simulation is fully periodic
 */
bool cfdemCloud::checkSimulationFullyPeriodic() {
  const polyBoundaryMesh& patches = mesh_.boundaryMesh();
  int nPatchesCyclic = 0;     // 周期边界数量
  int nPatchesNonCyclic = 0;  // 非周期边界数量
  for (const polyPatch& patch : patches) {
// 统计 nPatchesCyclic 和 nPatchesNonCyclic
#if defined(versionExt32)
    if (isA<cyclicPolyPatch>(patch)) {
      nPatchesCyclic += 1;
    } else if (!isA<processorPolyPatch>(patch)) {
      nPatchesNonCyclic += 1;
    }
#else
    if (isA<cyclicPolyPatch>(patch) || isA<cyclicAMIPolyPatch>(patch)) {
      nPatchesCyclic += 1;
    } else if (!isA<processorPolyPatch>(patch)) {
      nPatchesNonCyclic += 1;
    }
#endif
  }
  return nPatchesNonCyclic == 0;
}

tmp<volScalarField> cfdemCloud::ddtVoidFraction() const {
  if ("off" == ddtVoidFractionType()) {
    return tmp<volScalarField>(ddtVoidFraction_ * 0.);
  }
  return tmp<volScalarField>(ddtVoidFraction_ * 1.);
}

tmp<volScalarField> cfdemCloud::voidFractionNuEff(volScalarField& voidFraction) const {
  if (modelType() == "B" || modelType() == "Bfull" || modelType() == "none") {
#ifdef compre
    return tmp<volScalarField>(
        new volScalarField("viscousTerm", (turbulence_.mut() + turbulence_.mu() + turbulenceMultiphase_)));
#else
    return tmp<volScalarField>(
        new volScalarField("viscousTerm", (turbulence_.nut() + turbulence_.nu() + turbulenceMultiphase_)));
#endif
  } else {  // 如果是 model A, 则需要乘以空隙率
#ifdef compre
    return tmp<volScalarField>(new volScalarField(
        "viscousTerm", voidFraction * (turbulence_.mut() + turbulence_.mu() + turbulenceMultiphase_)));
#else
    return tmp<volScalarField>(new volScalarField(
        "viscousTerm", voidFraction * (turbulence_.nut() + turbulence_.nu() + turbulenceMultiphase_)));
#endif
  }
}

tmp<fvVectorMatrix> cfdemCloud::divVoidFractionTau(volVectorField& U, volScalarField& voidFraction) const {
  return -fvm::laplacian(voidFractionNuEff(voidFraction), U) -
         fvc::div(voidFractionNuEff(voidFraction) * dev2(fvc::grad(U)().T()));
}

// Foam::autoPtr<T> 中定义了 inline operator const T&() const;
const dataExchangeModel& cfdemCloud::dataExchangeM() const {
  return dataExchangeModel_;
}

const voidFractionModel& cfdemCloud::voidFractionM() const {
  return voidFractionModel_;
}

const locateModel& cfdemCloud::locateM() const {
  return locateModel_;
}

const averagingModel& cfdemCloud::averagingM() const {
  return averagingModel_;
}

const globalForce& cfdemCloud::globalF() const {
  return globalForce_;
}

dataExchangeModel& cfdemCloud::dataExchangeM() {
  return dataExchangeModel_();
}

voidFractionModel& cfdemCloud::voidFractionM() {
  return voidFractionModel_();
}

locateModel& cfdemCloud::locateM() {
  return locateModel_();
}

averagingModel& cfdemCloud::averagingM() {
  return averagingModel_();
}

globalForce& cfdemCloud::globalF() {
  return globalForce_();
}

#if defined(version24Dev)
const turbulenceModel& cfdemCloud::turbulence() const
#elif defined(version21) || defined(version16ext)
#ifdef compre
const compressible::turbulenceModel& cfdemCloud::turbulence() const
#else
const incompressible::turbulenceModel& cfdemCloud::turbulence() const
#endif
#elif defined(version15)
const incompressible::RASModel& cfdemCloud::turbulence() const
#endif
{
  return turbulence_;
}

}  // namespace Foam
