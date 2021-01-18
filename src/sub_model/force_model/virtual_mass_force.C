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

#include "./global_force.h"
#include "./virtual_mass_force.h"

namespace Foam {

cfdemDefineTypeName(virtualMassForce);

cfdemCreateNewFunctionAdder(forceModel, virtualMassForce);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
virtualMassForce::virtualMassForce(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      U_(cloud.globalF().U()),
      voidFraction_(cloud.globalF().voidFraction()),
      phi_(cloud.globalF().phi()) {
  createForceSubModels(subPropsDict_, kUnResolved);
}

virtualMassForce::~virtualMassForce() {}

void virtualMassForce::setForce() {
  Info << "Setting virtual mass force..." << endl;
  base::MPI_Barrier();
  UInterpolator_.reset(
      interpolation<Foam::vector>::New(subPropsDict_.lookupOrDefault("UInterpolationType", word("cellPointFace")), U_)
          .ptr());
  voidFractionInterpolator_.reset(
      interpolation<Foam::scalar>::New(
          subPropsDict_.lookupOrDefault("voidfractionInterpolationType", word("cellPoint")), voidFraction_)
          .ptr());
  DDtUInterpolator_.reset(
      interpolation<vector>::New(subPropsDict_.lookupOrDefault("DDtUInterpolationType", word("cellPointFace")),
                                 cloud_.globalF().ddtU())
          .ptr());
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // 只设置 fine and middle 颗粒的虚拟质量力
    if (cloud_.checkCoarseParticle(index)) {
      continue;
    }
    // 虚拟质量力
    Foam::vector virtualMassForce = Foam::vector::zero;
    setForceKernel(index, virtualMassForce);
    // setMiddleParticleForceKernel(virtualMassForce, index);
    forceSubModel_->partToArray(index, virtualMassForce, Foam::vector::zero, Foam::vector::zero, 0);
  }
  Info << "Setting virtual mass force - done" << endl;
  base::MPI_Barrier();
}

void virtualMassForce::setForceKernel(const int index, Foam::vector& virtualMassForce) {
  // 颗粒中心所在网格的索引
  int findCellID = cloud_.findCellIDs()[index];
  // 获取背景流体ddtU
  Foam::vector ddtU = getBackgroundDDtU(index, findCellID);
  // 获取背景流体速度
  Foam::vector Ufluid = getBackgroundUfluid(index, findCellID);
  // 获取背景网格空隙率
  double vf = getBackgroundVoidFraction(index, findCellID);
  // 对于每一个颗粒，只需要一个处理器计算，即颗粒中心所在的处理器
  if (findCellID >= 0) {
    double dt = cloud_.mesh().time().deltaT().value();     // 时间步长
    double radius = cloud_.getRadius(index);               // 颗粒半径
    double diameter = radius;                              // 颗粒直径
    double rho = cloud_.globalF().rhoField()[findCellID];  // 流体密度
    double Ac = 0.0;                                       // Ac 系数
    Foam::vector Up = cloud_.getVelocity(index);           // 颗粒速度
    Foam::vector prevUp = Foam::vector::zero;              // 上一个时间步颗粒速度
    Foam::vector ddtUp = Foam::vector::zero;               // 颗粒加速度
    Foam::vector Ur = Foam::vector::zero;                  // 相对速度
    Foam::vector ddtUr = Foam::vector::zero;               // 相对加速度
    // 计算颗粒加速度 ddtUp
    prevUp = cloud_.globalF().getPrevParticleVel(index);
    ddtUp = (Up - prevUp) / dt;
    // 计算 Ac 系数
    Ur = Ufluid - Up;
    ddtUr = ddtU - ddtUp;
    Ac = sqr(mag(Ur)) / (diameter * mag(ddtUr));
    // 计算虚拟质量力
    virtualMassForce = (2.0 / 3.0) * M_PI * (2.1 - 0.132 / (0.12 + Ac * Ac)) * pow(diameter, 3) * rho * ddtUr;
    if ("B" == cloud_.modelType()) {
      virtualMassForce /= vf;
    }
    if (forceSubModel_->verbose()) {
      Pout << "vf = " << vf << endl;
      Pout << "Ufluid = " << Ufluid << endl;
      Pout << "ddtU = " << ddtU << endl;
      Pout << "ddtUp = " << ddtUp << endl;
      Pout << "ddtUr = " << ddtUr << endl;
      Pout << "Ac = " << Ac << endl;
      Pout << "virtual mass force = " << virtualMassForce << endl;
    }
  }
}

//! \brief 计算颗粒 index 处的背景流体的ddtu
Foam::vector virtualMassForce::getBackgroundDDtU(const int index, const int findCellID) const {
  Foam::vector ddtU = Foam::vector::zero;
  // 对于 middle 颗粒，使用高斯核函数计算
  if (forceSubModel_->useGaussCoreFunctionRefined() && cloud_.checkMiddleParticle(index)) {
    ddtU = cloud_.globalF().getBackgroundDDtU(index);
  } else if (!cloud_.checkCoarseParticle(index) && findCellID >= 0) {
    if (forceSubModel_->interpolation()) {
      // 获取颗粒中心的坐标, 将颗粒中心所在网格的空隙率和流体速度插值到颗粒中心处
      Foam::vector pos = cloud_.getPosition(index);
      ddtU = DDtUInterpolator_().interpolate(pos, findCellID);
    } else {
      // 不使用插值模型，直接设置颗粒中心处的值为颗粒中心所在网格的值
      ddtU = cloud_.globalF().ddtU()[findCellID];
    }
  }
  return ddtU;
}

//! \brief 计算颗粒 index 处的背景流体速度
Foam::vector virtualMassForce::getBackgroundUfluid(const int index, const int findCellID) const {
  // 背景流体速度
  Foam::vector Ufluid = Foam::vector::zero;
  // 对于 middle 颗粒，使用高斯核函数计算
  if (forceSubModel_->useGaussCoreFunctionRefined() && cloud_.checkMiddleParticle(index)) {
    Ufluid = cloud_.globalF().getBackgroundUfluid(index);
  } else if (!cloud_.checkCoarseParticle(index) && findCellID >= 0) {
    if (forceSubModel_->interpolation()) {
      // 获取颗粒中心的坐标, 将颗粒中心所在网格的空隙率和流体速度插值到颗粒中心处
      Foam::vector pos = cloud_.getPosition(index);
      Ufluid = UInterpolator_().interpolate(pos, findCellID);
    } else {
      // 不使用插值模型，直接设置颗粒中心处的值为颗粒中心所在网格的值
      Ufluid = U_[findCellID];
    }
  }
  return Ufluid;
}

//! \brief 计算颗粒 index 处的背景流体空隙率
double virtualMassForce::getBackgroundVoidFraction(const int index, const int findCellID) const {
  double vf = 1.0;
  // 对于 middle 颗粒，使用高斯核函数计算
  if (forceSubModel_->useGaussCoreFunctionRefined() && cloud_.checkMiddleParticle(index)) {
    vf = cloud_.globalF().getBackgroundVoidFraction(index);
  } else if (!cloud_.checkCoarseParticle(index) && findCellID >= 0) {
    if (forceSubModel_->interpolation()) {
      // 获取颗粒中心的坐标, 将颗粒中心所在网格的空隙率和流体速度插值到颗粒中心处
      Foam::vector pos = cloud_.getPosition(index);
      vf = voidFractionInterpolator_().interpolate(pos, findCellID);
      // 确保插值后颗粒中心的空隙率有意义
      vf = vf > 1.0 ? 1.0 : vf;
      vf = vf < Foam::SMALL ? Foam::SMALL : vf;
    } else {
      // 不使用插值模型，直接设置颗粒中心处的值为颗粒中心所在网格的值
      vf = voidFraction_[findCellID];
    }
  }
  return vf;
}

void virtualMassForce::setMiddleParticleForceKernel(Foam::vector& virtualMassForce, const int index) {
  int findCellID = cloud_.findCellIDs()[index];
  if (findCellID >= 0) {
    std::unordered_set<int> set;                           // 颗粒扩展网格的集合
    double dt = cloud_.mesh().time().deltaT().value();     // 时间步长
    double radius = cloud_.getRadius(index);               // 颗粒半径
    double diameter = radius;                              // 颗粒直径
    double rho = cloud_.globalF().rhoField()[findCellID];  // 流体密度
    double vf = 0.0;                                       // 流体空隙率
    double Ac = 0.0;                                       // Ac 系数
    Foam::vector particlePos = cloud_.getPosition(index);  // 颗粒中心坐标
    Foam::vector Up = cloud_.getVelocity(index);           // 颗粒速度
    Foam::vector prevUp = Foam::vector::zero;              // 上一个时间步颗粒速度
    Foam::vector UFluid = Foam::vector::zero;              // 背景流体速度
    Foam::vector ddtU = Foam::vector::zero;                // 背景流体加速度
    Foam::vector ddtUp = Foam::vector::zero;               // 颗粒加速度
    Foam::vector Ur = Foam::vector::zero;                  // 相对速度
    Foam::vector ddtUr = Foam::vector::zero;               // 相对加速度
    // 获取当前颗粒的 Expanded Cell 集合
    cloud_.voidFractionM().buildExpandedCellSet(set, findCellID, particlePos, radius, 6);
    // 计算背景流体速度与背景流体加速度
    double sumCore = 0.0;
    double core = 0.0;
    double cellV = 0.0;
    double sumPV = 0.0;
    double sumCV = 0.0;
    Foam::vector cellPos = Foam::vector::zero;
    for (int cellID : set) {
      if (cellID >= 0) {  // cell found
        cellPos = cloud_.mesh().C()[cellID];
        cellV = cloud_.mesh().V()[cellID];
        // 计算高斯核
        core = globalForce::GaussCore(particlePos, cellPos, radius, 6);
        // 计算累计速度
        UFluid += voidFraction_[cellID] * U_[cellID] * core * cellV;
        // 计算累计加速度
        ddtU += voidFraction_[cellID] * cloud_.globalF().ddtU()[cellID] * core * cellV;
        // 计算累加因数
        sumCore += voidFraction_[cellID] * core * cellV;
        // 计算累加流体体积
        sumPV += voidFraction_[cellID] * cellV;
        // 计算累加网格体积
        sumCV += cellV;
      }
    }
    UFluid /= sumCore;
    ddtU /= sumCore;
    vf = sumPV / sumCV;
    // 计算颗粒加速度 ddtUp
    prevUp = cloud_.globalF().getPrevParticleVel(index);
    ddtUp = (Up - prevUp) / dt;
    // 计算 Ac 系数
    Ur = UFluid - Up;
    ddtUr = ddtU - ddtUp;
    Ac = sqr(mag(Ur)) / (diameter * mag(ddtUr));
    // 计算虚拟质量力
    virtualMassForce = (2.0 / 3.0) * M_PI * (2.1 - 0.132 / (0.12 + Ac * Ac)) * pow(diameter, 3) * rho * ddtUr;
    if ("B" == cloud_.modelType()) {
      virtualMassForce /= vf;
    }
    if (forceSubModel_->verbose()) {
      Pout << "vf = " << vf << endl;
      Pout << "ddtU = " << ddtU << endl;
      Pout << "ddtUp = " << ddtUp << endl;
      Pout << "ddtUr = " << ddtUr << endl;
      Pout << "Ac = " << Ac << endl;
      Pout << "virtual mass force = " << virtualMassForce << endl;
    }
  }
}

}  // namespace Foam
