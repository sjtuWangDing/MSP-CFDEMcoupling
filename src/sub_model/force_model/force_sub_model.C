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
  This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
  and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
  forceSubModel
\*---------------------------------------------------------------------------*/

#include "./force_model.h"
#include "./force_sub_model.h"

namespace Foam {

const int forceSubModel::Switches::kNum = 7;

const char* forceSubModel::Switches::kNameList[] = {"treatForceExplicitInMomEquation",
                                                    "treatForceBothCFDAndDEM",
                                                    "treatDEMForceImplicit",
                                                    "verbose",
                                                    "interpolation",
                                                    "scalarViscosity",
                                                    "useGaussCoreFunctionRefined"};

//! \brief Constructor
forceSubModel::forceSubModel(cfdemCloud& cloud, forceModel& forceModel, const dictionary& subPropsDict)
    : cloud_(cloud),
      forceModel_(forceModel),
      subPropsDict_(subPropsDict),
      switches_(),
      rho_(cloud.globalF().rhoField()),
      nu_(IOobject("scalarViscosity", cloud.mesh().time().timeName(), cloud.mesh(), IOobject::NO_READ,
                   IOobject::NO_WRITE),
          cloud.mesh(), dimensionedScalar("nu0", dimensionSet(0, 2, -1, 0, 0), 1.0)) {}

//! \brief Destructor
forceSubModel::~forceSubModel() {}

/*!
 * \param index                  <[in] 颗粒索引
 * \param dragTot                <[in] 索引为 index 的颗粒受到的总阻力
 * \param dragEx                 <[in] 索引为 index 的颗粒受到的显式阻力
 * \param Ufluid = vector::zero  <[in] 索引为 index 的颗粒中心处流体速度(可以指定是否使用插值模型计算)
 * \param scalar Cd = 0          <[in] 颗粒阻力系数
 */
void forceSubModel::partToArray(const int& index, const Foam::vector& dragTot, const Foam::vector& dragEx,
                                const Foam::vector& Ufluid, scalar Cd) const {
  if (treatForceBothCFDAndDEM()) {  // CFD 与 DEM 求解器都考虑耦合力
    if (treatForceExplicitInMomEquation()) {
#pragma unroll
      // 耦合力为显式力
      for (int j = 0; j < 3; ++j) {
        // 将耦合力累加到 expFoces_
        cloud_.expForces()[index][j] += dragTot[j];
      }
    } else {  // 耦合力为隐式力
#pragma unroll
      for (int j = 0; j < 3; ++j) {
        // 将耦合力累加到 impForces_ 和 expFoces_
        cloud_.impForces()[index][j] += dragTot[j] - dragEx[j];  // 隐式力 = dragTot[j] - dragEx[j]
        cloud_.expForces()[index][j] += dragEx[j];
      }
    }
  }
  if (treatDEMForceImplicit()) {
#pragma unroll
    // 颗粒中心处的流体速度和阻力系数都被传递到 DEM 中，从而在每个 DEM
    // 时间步中，使用阻力系数和流体速度，与当前颗粒速度一起计算颗粒受到的阻力
    for (int j = 0; j < 3; j++) {
      cloud_.fluidVel()[index][j] = Ufluid[j];
    }
    cloud_.cds()[index] = Cd;
  } else {
#pragma unroll
    // 直接将总阻力传递给 DEMForces_[index]
    // usually used for ArchimedesIB and ShirgaonkarIB force model
    for (int j = 0; j < 3; j++) {
      // 如果使用 Junke Guo BBO Equation，则需要 3 倍的浮力
      if (forceModel_.typeName() == "Archimedes" && cloud_.useGuoBBOEquation()) {
        cloud_.DEMForces()[index][j] += 3.0 * dragTot[j];
      } else {
        cloud_.DEMForces()[index][j] += dragTot[j];
      }
    }
  }
}

/*!
 * \param index                  <[in] 颗粒索引
 * \param torque = vector::zero  <[in] 索引为 index 的颗粒受到的力矩
 */
void forceSubModel::addTorque(int index, const Foam::vector& torque /* = Foam::vector::zero */) const {
  for (int i = 0; i < 3; ++i) {
    cloud_.DEMTorques()[index][i] += torque[i];
  }
}

//! \brief read switches from force model dictionary
void forceSubModel::readSwitches() {
  // 遍历每一个 switch
  for (int i = 0; i < Switches::kNum; ++i) {
    Info << "Foam::forceSubModel::readSwitches(): looking for " << Switches::kNameList[i] << "...";
    if (subPropsDict_.found(Switches::kNameList[i])) {
      bool res = Switch(subPropsDict_.lookup(Switches::kNameList[i]));
      Info << " found in dict, using " << (res ? "true" : "false") << endl;
      // 设置当前的 switch
      if (res) {
        switches_.setTrue(static_cast<ESwitch>(i));
      }
    } else {
      // 如果在 dict 中没有找到，则使用默认值 false
      Info << " not found in dict, using false" << endl;
    }
  }
  // try to read scalarViscosity
  if (switches_.isTrue(kScalarViscosity)) {
    dimensionedScalar nu0_("nu", dimensionSet(0, 2, -1, 0, 0), subPropsDict_.lookup("nu"));
    nu_ = volScalarField(IOobject("scalarViscosity", cloud_.mesh().time().timeName(), cloud_.mesh(), IOobject::NO_READ,
                                  IOobject::NO_WRITE),
                         cloud_.mesh(), nu0_);
  }
  // TODO: 检查 switch 是否存在冲突
}

/*!
 * \brief check switches for different type of force
 * \param forceType force type, Eg: kUnResolved
 */
void forceSubModel::checkSwitches(EForceType forceType) const {
  switch (forceType) {
    case kUnResolved:
      break;
    case kSemiResolved:
      break;
    case kResolved:
      // for resolved force, kTreatForceBothCFDAndDEM == false and kTreatDEMForceImplicit == false
      if (switches_.isTrue(kTreatForceBothCFDAndDEM) || switches_.isTrue(kTreatDEMForceImplicit)) {
        FatalError << "Error: illegal force switches for kResolved force type" << abort(FatalError);
      }
      break;
    case kMix:
      break;
    default:
      FatalError << "Error: illegal force type: " << forceType << abort(FatalError);
  }
}

/*!
 * \brief 计算 IB drag，用于计算 ShirgaonkarIBModel and mixShirgaonkarIBModel
 *        dimensionSet(1, -2, -2, 0, 0)
 * \param U 速度场
 * \param p 压力场
 * \return IB drag
 */
volVectorField forceSubModel::IBDrag(const volVectorField& U, const volScalarField& p) const {
#ifdef compre
  return muField() * fvc::laplacian(U) - fvc::grad(p);
#else
  return rhoField() * (nuField() * fvc::laplacian(U) - fvc::grad(p));
#endif
}

volVectorField forceSubModel::divTauField(const volVectorField& U) const {
#ifdef compre
  const volScalarField& mu_ = muField();
  return -fvc::laplacian(mu_, U) - fvc::div(mu_ * dev(fvc::grad(U)().T()));
#else
  const volScalarField& nu_ = nuField();
  const volScalarField& rho_ = rhoField();
  return -fvc::laplacian(nu_ * rho_, U) - fvc::div(nu_ * rho_ * dev(fvc::grad(U)().T()));
#endif
}

}  // namespace Foam
