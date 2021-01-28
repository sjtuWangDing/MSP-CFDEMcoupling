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
\*---------------------------------------------------------------------------*/

#include "./Basset_force.h"
#include "./global_force.h"
#include "./virtual_mass_force.h"
#include "cfdem_tools/cfdem_tools.h"

namespace Foam {

class defaultField {
 public:
  static const volScalarField& getVoidFraction(const fvMesh& mesh) {
    if (!voidFractionSp_) {
      voidFractionSp_ = std::make_shared<volScalarField>(
          IOobject("voidFractionDefault", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh,
          dimensionedScalar("1", dimensionSet(0, 0, 0, 0, 0), 1.0));
    }
    return *voidFractionSp_;
  }
  static const volScalarField& getVolumeFraction(const fvMesh& mesh) {
    if (!volumeFractionSp_) {
      volumeFractionSp_ = std::make_shared<volScalarField>(
          IOobject("volumeFractionDefault", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh,
          dimensionedScalar("1", dimensionSet(0, 0, 0, 0, 0), 1.0));
    }
    return *volumeFractionSp_;
  }

 private:
  static std::shared_ptr<volScalarField> voidFractionSp_;
  static std::shared_ptr<volScalarField> volumeFractionSp_;
};

std::shared_ptr<volScalarField> defaultField::voidFractionSp_(nullptr);
std::shared_ptr<volScalarField> defaultField::volumeFractionSp_(nullptr);

cfdemDefineTypeName(globalForce);

cfdemDefineNewFunctionMap(globalForce);

cfdemDefineConstructNewFunctionMap(globalForce);

cfdemDefineDestroyNewFunctionMap(globalForce);

cfdmeDefineBaseTypeNew(autoPtr, globalForce, (cfdemCloud & cloud, const dictionary& dict), dict, (cloud));

cfdemCreateNewFunctionAdder(globalForce, globalForce);

globalForce::globalForce(cfdemCloud& cloud)
    : cloud_(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      verbose_(subPropsDict_.lookupOrDefault<bool>("verbose", false)),
      GaussCoreEff_(subPropsDict_.lookupOrDefault<double>("GaussCoreEff", 3.0)),
      ddtU_(IOobject("ddtU", cloud.mesh().time().timeName(), cloud.mesh(), IOobject::READ_IF_PRESENT,
                     IOobject::AUTO_WRITE),
            cloud.mesh(),
            dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector(0, 0, 0))),  // [ddtU] == [m / s^2])
      impParticleForce_(IOobject("impParticleForce", cloud.mesh().time().timeName(), cloud.mesh(),
                                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                        cloud.mesh(), dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0),
                                                        vector(0, 0, 0))),  // [N] == [kg * m / s^2]
      expParticleForce_(IOobject("expParticleForce", cloud.mesh().time().timeName(), cloud.mesh(),
                                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                        cloud.mesh(), dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0),
                                                        vector(0, 0, 0))),  // [N] == [kg * m / s^2]
      gravityFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("gravityFieldName", "g").c_str()),
      densityFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("densityFieldName", "rho").c_str()),
      velFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("velFieldName", "U").c_str()),
      pressureFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("pressureFieldName", "p").c_str()),
      phiFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("phiFieldName", "phi").c_str()),
      voidFractionFieldName_(
          subPropsDict_.lookupOrDefault<Foam::word>("voidFractionFieldName", "voidFraction").c_str()),
      volumeFractionFieldName_(
          subPropsDict_.lookupOrDefault<Foam::word>("volumeFractionFieldName", "volumeFraction").c_str()),
#if defined(version21)
      g_(cloud.mesh().lookupObject<uniformDimensionedVectorField>(gravityFieldName_)),
#elif defined(version16ext) || defined(version15)
      g_(dimensionedVector(
             cloud.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(environmentalProperties))
             .value()),
#endif
      rho_(cloud.mesh().lookupObject<volScalarField>(densityFieldName_)),
      U_(cloud.mesh().lookupObject<volVectorField>(velFieldName_)),
      p_(cloud.mesh().lookupObject<volScalarField>(pressureFieldName_)),
      phi_(cloud.mesh().lookupObject<surfaceScalarField>(phiFieldName_)),
      voidFraction_(subPropsDict_.found("voidFractionFieldName")
                        ? cloud.mesh().lookupObject<volScalarField>(voidFractionFieldName_)
                        : defaultField::getVoidFraction(cloud.mesh())),
      volumeFraction_(subPropsDict_.found("volumeFractionFieldName")
                          ? cloud.mesh().lookupObject<volScalarField>(volumeFractionFieldName_)
                          : defaultField::getVolumeFraction(cloud.mesh())) {
}

globalForce::~globalForce() {}

//! \brief 每一次耦合中，在 set force 前执行
void globalForce::initBeforeSetForce() {
  bool isUsedVirtualMassForce = cfdemTools::isUsedForceModel(cloud_, virtualMassForce::cTypeName());
  bool isUsedBassetForce = cfdemTools::isUsedForceModel(cloud_, BassetForce::cTypeName());
  // 如果使用了 BassetForce or virtualMassForce，则需要计算 ddtU
  if (isUsedVirtualMassForce || isUsedBassetForce) {
    // 计算 ddtU field
    ddtU_ = fvc::ddt(U_) + fvc::div(phi_, U_);
  }
  base::MPI_Info("globalForce: initBeforeSetForce - done", verbose_);
}

//! \brief 每一次耦合中，在 set force 后执行
void globalForce::endAfterSetForce() {
  bool isUsedVirtualMassForce = cfdemTools::isUsedForceModel(cloud_, virtualMassForce::cTypeName());
  bool isUsedBassetForce = cfdemTools::isUsedForceModel(cloud_, BassetForce::cTypeName());
  // 如果使用了 BassetForce or virtualMassForce，则需要保存颗粒速度以及更新 ddtUrHistoryMap_
  if (isUsedVirtualMassForce || isUsedBassetForce) {
    for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
      // (1) 更新颗粒速度，所有处理器都必须更新，当颗粒从一个计算域运动到另一个计算域，新的计算域必须有颗粒的 prevU
      updatePrevParticleVelMap(index);
      // (2) 更新 ddtUrHistoryMap_，所有处理器都必须更新
      if (isUsedBassetForce) {
        updateDDtUrHistoryMap(index);
      }
    }  // end of loop particles
  }
  base::MPI_Info("globalForce: endAfterSetForce - done", verbose_);
}

//! \brief 更新颗粒速度
void globalForce::updatePrevParticleVelMap(const int index) {
  if (prevParticleVelMap_.end() == prevParticleVelMap_.find(index)) {
    prevParticleVelMap_.insert(std::make_pair(index, cloud_.getVelocity(index)));
  } else {
    prevParticleVelMap_[index] = cloud_.getVelocity(index);
  }
}

//! \brief 更新 ddtUrHistoryMap
void globalForce::updateDDtUrHistoryMap(const int index) {
  int procId = base::procId();                         // 处理器编号
  int numProc = base::numProc();                       // 处理器数量
  int rootProc = cloud_.particleRootProcIDs()[index];  // 主节点编号，即颗粒所在的处理器编号
  // 只有一个节点直接跳过
  if (1 == numProc) {
    return;
  }
  std::vector<Foam::vector>& ddtUrHistoryVec = getDDtUrHistory(index);  // ddtUrHistory
  std::vector<double> ddtUr(3, 0.0);                                    // data buffer
  // 颗粒中心所在的主节点作为广播节点填充 data buffer
  if (rootProc == procId) {
    CHECK(!ddtUrHistoryVec.empty()) << __func__ << ": ddtUrHistoryVec of particle " << index
                                    << " in the root proc is empty";
    int size = ddtUrHistoryVec.size();
    ddtUr[0] = ddtUrHistoryVec[size - 1][0];
    ddtUr[1] = ddtUrHistoryVec[size - 1][1];
    ddtUr[2] = ddtUrHistoryVec[size - 1][2];
  }
  // rootProc 节点广播数据
  MPI_Bcast(ddtUr.data(), ddtUr.size(), MPI_DOUBLE, rootProc, MPI_COMM_WORLD);
  // 所有的非主节点保存数据
  if (rootProc != procId) {
    ddtUrHistoryVec.push_back(Foam::vector(ddtUr[0], ddtUr[1], ddtUr[2]));
  }
  base::MPI_Barrier();
}

//! \brief 获取上一个耦合时间步中颗粒速度
Foam::vector globalForce::getPrevParticleVel(const int index) const {
  auto iter = prevParticleVelMap_.find(index);
  if (prevParticleVelMap_.end() != iter) {
    return iter->second;
  } else if (cloud_.dataExchangeM().isFirstCouplingStep()) {
    // 如果是第一个耦合时间步，则使用初始颗粒速度
    return cloud_.getInitVelocity(index);
  }
  return Foam::vector::zero;
}

//! \brief 获取颗粒的历史 ddtUr
std::vector<Foam::vector>& globalForce::getDDtUrHistory(const int index) {
  auto iter = ddtUrHistoryMap_.find(index);
  if (ddtUrHistoryMap_.end() != iter) {
    // find particle of index's history ddtUr
    return iter->second;
  }
  // not find
  ddtUrHistoryMap_.insert(std::make_pair(index, std::vector<Foam::vector>()));
  return ddtUrHistoryMap_[index];
}

}  // namespace Foam
