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
#include "sub_model/averaging_model/averaging_model.h"
#include "sub_model/data_exchange_model/data_exchange_model.h"
#include "sub_model/force_model/force_model.h"
#include "sub_model/liggghts_command_model/liggghts_command_model.h"
#include "sub_model/locate_model/locate_model.h"
#include "sub_model/void_fraction_model/void_fraction_model.h"

namespace Foam {

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
      dataExchangeModel_(dataExchangeModel::New(*this, couplingPropertiesDict_)),
      voidFractionModel_(voidFractionModel::New(*this, couplingPropertiesDict_)),
      locateModel_(locateModel::New(*this, couplingPropertiesDict_)),
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
  // check periodic
  if (checkPeriodicCells() != checkSimulationFullyPeriodic()) {
    FatalError << "checkSimulationFullyPeriodic(): " << (checkSimulationFullyPeriodic() ? true : false)
               << ", but from dictionary read checkPeriodicCells: " << (checkPeriodicCells() ? true : false)
               << abort(FatalError);
  }
}

cfdemCloud::~cfdemCloud() {}

//! \brief 重新分配内存
void cfdemCloud::reallocate() {
  // if (numberOfParticlesChanged()) {
  //   int number = numberOfParticles();
  //   // allocate memory of data exchanged with liggghts
  //   dataExchangeM().realloc(parCloud_.radii(), base::makeShape1(number),
  //   parCloud_.radiiPtr(), 0.0);
  //   dataExchangeM().realloc(parCloud_.positions(), base::makeShape2(number,
  //   3), parCloud_.positionsPtr(), 0.0);
  //   dataExchangeM().realloc(parCloud_.velocities(), base::makeShape2(number,
  //   3), parCloud_.velocitiesPtr(), 0.0);
  //   // allocate memory of data not exchanged with liggghts
  //   parCloud_.particleOverMeshNumber() =
  //   std::move(base::CITensor1(base::makeShape1(number), 0));
  //   parCloud_.findCellIDs() =
  //   std::move(base::CITensor1(base::makeShape1(number), -1));
  //   parCloud_.dimensionRatios() =
  //   std::move(base::CDTensor1(base::makeShape1(number), 0.0));
  // }
}

/*!
 * \brief 更新函数
 * \note used for cfdemSolverPiso
 * \param VoidF  <[in, out] 小颗粒空隙率场
 * \param Us     <[in, out] 局部平均小颗粒速度场
 * \param U      <[in] 流体速度场
 */
void cfdemCloud::evolve(volScalarField& VoidF, volVectorField& Us, volVectorField& U) {
  // Info << "\nFoam::cfdemCloud::evolve(), used for cfdemSolverPiso......\n" <<
  // endl;
  // if (!writeTimePassed_ && mesh_.time().outputTime()) {
  //   writeTimePassed_ = true;
  // }
  // if (dataExchangeM().doCoupleNow()) {
  //   Info << "evolve coupling..." << endl;
  //   // couple() 函数执行 liggghts 脚本，并获取新的颗粒数量
  //   parCloud_.setNumberOfParticles(dataExchangeM().couple());
  //   Info << "get number of particles: " << parCloud_.numberOfParticles()<< "
  //   at coupling step: "
  //     << dataExchangeM().couplingStep() << endl;
  //   // 重置局部平均颗粒速度
  //   averagingM().resetUs();
  //   // 重置颗粒速度影响因数场
  //   averagingM().resetUsWeightField();
  //   // 重置小颗粒空隙率场
  //   voidFractionM().resetVoidFraction();
  //   // 重置隐式力场
  //   // 重置显式力场
  //   // 重置动量交换场
  // }
  // Info << "Foam::cfdemCloud::evolve() - done\n" << endl;
  // return true;
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

}  // namespace Foam
