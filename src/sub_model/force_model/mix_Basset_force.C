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
#include "./mix_Basset_force.h"

namespace Foam {

cfdemDefineTypeName(mixBassetForce);

cfdemCreateNewFunctionAdder(forceModel, mixBassetForce);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
mixBassetForce::mixBassetForce(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      U_(cloud.globalF().U()),
      voidFraction_(cloud.globalF().voidFraction()),
      phi_(cloud.globalF().phi()) {
  createForceSubModels(subPropsDict_, kUnResolved);
}

mixBassetForce::~mixBassetForce() {}

void mixBassetForce::setForce() {
  Info << "Setting mix Basset force..." << endl;
  base::MPI_Barrier();
  if (cloud_.numberOfParticlesChanged()) {
    FatalError << "Currently mixBassetForce model not allow number of particle changed\n" << abort(FatalError);
  }
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // 只设置 fine and middle 颗粒的 Basset force
    if (cloud_.checkCoarseParticle(index)) {
      continue;
    }
    // Basset force
    Foam::vector BstForce = Foam::vector::zero;
    setForceKernel(index, BstForce);
    forceSubModel_->partToArray(index, BstForce, Foam::vector::zero, Foam::vector::zero, 0.0);
  }
  Info << "Setting mix Basset force - done" << endl;
  base::MPI_Barrier();
}

void mixBassetForce::setForceKernel(const int index, Foam::vector& BstForce) {
  // 颗粒中心所在网格的索引
  int findCellID = cloud_.findCellIDs()[index];
  // 对于每一个颗粒，只需要一个处理器计算，即颗粒中心所在的处理器
  if (findCellID >= 0) {
    double dt = cloud_.mesh().time().deltaT().value();                 // 时间步长
    double radius = cloud_.getRadius(index);                           // 颗粒半径
    double diameter = 2 * radius;                                      // 颗粒直径
    double rho = cloud_.globalF().rhoField()[findCellID];              // 流体密度
    double nu = forceSubModel_->nuField()[findCellID];                 // 流体动力粘度
    double B = 1.5 * sqr(diameter) * sqrt(M_PI * rho * rho * nu);      // Basset force 系数
    double vf = getBackgroundVoidFraction(index, findCellID);          // 背景网格空隙率
    Foam::vector ddtU = getBackgroundDDtU(index, findCellID);          // 背景流体ddtU
    Foam::vector Up = cloud_.getVelocity(index);                       // 颗粒速度
    Foam::vector prevUp = cloud_.globalF().getPrevParticleVel(index);  // 上一个时间步颗粒速度
    Foam::vector ddtUp = Foam::vector::zero;                           // 颗粒加速度
    Foam::vector ddtUr = Foam::vector::zero;                           // 相对加速度
    Foam::vector ddtUrHistory = Foam::vector::zero;                    // 累计相对加速度
    // 计算颗粒加速度 ddtUp
    ddtUp = (Up - prevUp) / dt;
    // 计算当前时间步的 ddtUr
    ddtUr = ddtU - ddtUp;
    // 计算 ddtUr 的历史累计和
    std::vector<Foam::vector>& ddtUrVec = cloud_.globalF().getDDtUrHistory(index);
    // 检查耦合时间步
    int size = ddtUrVec.size();
    CHECK_EQ(cloud_.dataExchangeM().couplingStep(), size) << __func__ << ": getDDtUrHistory record is error";
    // 计算 ddtUrHistory
    if (1 == size) {
      ddtUrHistory = 2 * ddtUrVec[0] / sqrt(dt);
    } else if (2 == size) {
      ddtUrHistory = ddtUrVec[0] / sqrt(2 * dt) + ddtUrVec[1] / sqrt(dt);
    } else if (2 < size) {
      ddtUrHistory = ddtUrVec[0] / sqrt(size * dt) + ddtUrVec[size - 1] / sqrt(1 * dt);
      for (int i = 1; i < size - 1; ++i) {
        ddtUrHistory += 2 * ddtUrVec[i] / sqrt((size - i) * dt);
      }
    }  // no need to calcualte ddtUrHistory when 0 == size
#if 1
    // 计算 Basset forece
    // TODO:
    // BstForce = 0.5 * B * dt * ddtUrHistory + 2 * B * ddtUr * sqrt(dt) +
    //            B * initUr / sqrt((cloud_.dataExchangeM().couplingStep() + 1) * dt);
    BstForce = 0.5 * B * dt * ddtUrHistory + 2 * B * ddtUr * sqrt(dt);
#else
    BstForce = 0.5 * B * dt * ddtUrHistory;
#endif
    // insert ddtUr to ddtUrVec
    ddtUrVec.push_back(ddtUr);
    if ("B" == cloud_.modelType()) {
      BstForce /= vf;
    }
    if (forceSubModel_->verbose()) {
      Pout << "ddtU = " << ddtU << endl;
      Pout << "Up = " << Up << endl;
      Pout << "prevUp = " << prevUp << endl;
      Pout << "ddtUp = " << ddtUp << endl;
      Pout << "ddtUr = " << ddtUr << endl;
      Pout << "mix Basset force = " << BstForce << endl;
    }
  }
}

//! \brief 计算颗粒 index 处的背景流体的ddtu
Foam::vector mixBassetForce::getBackgroundDDtU(const int index, const int findCellID) const {
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
Foam::vector mixBassetForce::getBackgroundUfluid(const int index, const int findCellID) const {
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
double mixBassetForce::getBackgroundVoidFraction(const int index, const int findCellID) const {
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

}  // namespace Foam
