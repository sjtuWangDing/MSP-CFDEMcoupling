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

namespace Foam {

cfdemDefineTypeName(BassetForce);

cfdemCreateNewFunctionAdder(forceModel, BassetForce);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
BassetForce::BassetForce(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      U_(cloud.globalF().U()),
      voidFraction_(cloud.globalF().voidFraction()),
      phi_(cloud.globalF().phi()) {
  createForceSubModels(subPropsDict_, kUnResolved);
}

BassetForce::~BassetForce() {}

//! \brief 更新颗粒速度，所有处理器都必须更新，否则当颗粒从一个计算域运动到另一个计算域，新的计算域必须有颗粒的 prevU
void BassetForce::updatePrevUp(const int index) {
  if (prevUpMap_.end() == prevUpMap_.find(index)) {
    prevUpMap_.insert(std::make_pair(index, cloud_.getVelocity(index)));
  } else {
    prevUpMap_[index] = cloud_.getVelocity(index);
  }
}

//! \brief 更新 DDtUrHistory
void BassetForce::updateDDtUrHistory(const int index) {
  Pout << "DDtUrHistoryMap_: " << DDtUrHistoryMap_[0].size() << endl;
}

void BassetForce::setForce() {
  Info << "Setting Basset force..." << endl;
  base::MPI_Barrier();
  if (cloud_.numberOfParticlesChanged()) {
    FatalError << "Currently BassetForce model not allow number of particle changed\n" << abort(FatalError);
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
  Info << "Setting Basset force - done" << endl;
  base::MPI_Barrier();
}

void BassetForce::setForceKernel(const int index, Foam::vector& BstForce) {
  // 颗粒中心所在网格的索引
  int findCellID = cloud_.findCellIDs()[index];
  // 获取背景流体ddtU
  Foam::vector ddtU = getBackgroundDDtU(index, findCellID);
  // 获取背景网格空隙率
  double vf = getBackgroundVoidFraction(index, findCellID);
  // 对于每一个颗粒，只需要一个处理器计算，即颗粒中心所在的处理器
  if (findCellID >= 0) {
    double dt = cloud_.mesh().time().deltaT().value();             // 时间步长
    double radius = cloud_.getRadius(index);                       // 颗粒半径
    double diameter = radius;                                      // 颗粒直径
    double rho = cloud_.globalF().rhoField()[findCellID];          // 流体密度
    double nu = forceSubModel_->nuField()[findCellID];             // 流体动力粘度
    double B = 1.5 * sqr(diameter) * sqrt(M_PI * rho * rho * nu);  // Basset force 系数
    Foam::vector Up = cloud_.getVelocity(index);                   // 颗粒速度
    Foam::vector initUr = cloud_.globalF().getInitUr(index);       // 颗粒初始时刻的相对速度
    Foam::vector prevUp = Foam::vector::zero;                      // 上一个时间步颗粒速度
    Foam::vector ddtUp = Foam::vector::zero;                       // 颗粒加速度
    Foam::vector ddtUr = Foam::vector::zero;                       // 相对加速度
    Foam::vector ddtUrHistory = Foam::vector::zero;                // 累计相对加速度
    // 计算颗粒加速度 ddtUp
    prevUp = cloud_.globalF().getPrevParticleVel(index);
    ddtUp = (Up - prevUp) / dt;
    // 计算当前时间步的 ddtUr
    ddtUr = ddtU - ddtUp;
    // 计算 ddtUr 的历史累计和
    std::vector<Foam::vector>& ddtUrVec = cloud_.globalF().getDDtUrHistory(index);
    int size = ddtUrVec.size();
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
    // insert ddtUr to ddtUrVec
    ddtUrVec.push_back(ddtUr);
    // 计算 Basset forece
    BstForce = 0.5 * B * dt * ddtUrHistory + 2 * B * ddtUr * sqrt(dt) +
               B * initUr / sqrt((cloud_.dataExchangeM().couplingStep() + 1) * dt);
    if ("B" == cloud_.modelType()) {
      BstForce /= vf;
    }
    if (forceSubModel_->verbose()) {
      Pout << "Basset force = " << BstForce << endl;
    }
  }
}

//! \brief 计算颗粒 index 处的背景流体的ddtu
Foam::vector BassetForce::getBackgroundDDtU(const int index, const int findCellID) const {
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
Foam::vector BassetForce::getBackgroundUfluid(const int index, const int findCellID) const {
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
double BassetForce::getBackgroundVoidFraction(const int index, const int findCellID) const {
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

#if 0
void BassetForce::setForce() {
  Info << "Setting Basset force..." << endl;
  if (cloud_.numberOfParticlesChanged()) {
    FatalError << "Currently BassetForce model not allow number of particle changed\n" << abort(FatalError);
  }
  const volScalarField& rhoField = forceSubModel_->rhoField();
  const volScalarField& nuField = forceSubModel_->nuField();
  // 计算 ddtU
  volVectorField ddtUField = fvc::ddt(U_) + fvc::div(phi_, U_);
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // Basset force
    Foam::vector BstForce = Foam::vector::zero;
    if (cloud_.checkMiddleParticle(index)) {
      setMiddleParticleForceKernel(BstForce, index, rhoField, nuField, ddtUField);
    }
    // 更新 Data
    updatePrevUp(index);
    updateDDtUrHistory(index);
    forceSubModel_->partToArray(index, BstForce, Foam::vector::zero, Foam::vector::zero, 0.0);
  }
  base::MPI_Barrier();
  Info << "Setting Basset force - done" << endl;
}

void BassetForce::setMiddleParticleForceKernel(Foam::vector& BstForce, const int index, const volScalarField& rhoField,
                                               const volScalarField& nuField, const volVectorField& ddtUField) {
  int findCellID = cloud_.findCellIDs()[index];
  if (findCellID >= 0) {
    std::unordered_set<int> set;                                   // 颗粒扩展网格的集合
    double dt = cloud_.mesh().time().deltaT().value();             // 时间步长
    double radius = cloud_.getRadius(index);                       // 颗粒半径
    double diameter = radius;                                      // 颗粒直径
    double rho = rhoField[findCellID];                             // 流体密度
    double nu = nuField[findCellID];                               // 流体动力粘度
    double vf = 0.0;                                               // 流体空隙率
    double B = 1.5 * sqr(diameter) * sqrt(M_PI * rho * rho * nu);  // Basset force 系数
    Foam::vector particlePos = cloud_.getPosition(index);          // 颗粒中心坐标
    Foam::vector Up = cloud_.getVelocity(index);                   // 颗粒速度
    Foam::vector prevUp = Foam::vector::zero;                      // 上一个时间步颗粒速度
    Foam::vector Ufluid = Foam::vector::zero;                      // 背景流体速度
    Foam::vector ddtU = Foam::vector::zero;                        // 背景流体加速度
    Foam::vector ddtUp = Foam::vector::zero;                       // 颗粒加速度
    Foam::vector ddtUr = Foam::vector::zero;                       // 相对加速度
    Foam::vector DDtUrHistory = Foam::vector::zero;                // 累计相对加速度

    // 获取当前颗粒的 Expanded Cell 集合
    cloud_.voidFractionM().buildExpandedCellSet(set, findCellID, particlePos, radius, 6);
    // 计算背景流体加速度
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
        Ufluid += voidFraction_[cellID] * U_[cellID] * core * cellV;
        // 计算累计加速度
        ddtU += voidFraction_[cellID] * ddtUField[cellID] * core * cellV;
        // 计算累加因数
        sumCore += voidFraction_[cellID] * core * cellV;
        // 计算累加流体体积
        sumPV += voidFraction_[cellID] * cellV;
        // 计算累加网格体积
        sumCV += cellV;
      }
    }
    // 计算平均流体速度
    Ufluid /= sumCore;
    if (cloud_.dataExchangeM().isFirstCouplingStep()) {
      initUr_ = (Ufluid - Up);
    }
    // 计算背景流体加速度
    ddtU /= sumCore;
    // 计算空隙率
    vf = sumPV / sumCV;
    // 计算颗粒加速度 ddtUp
    if (prevUpMap_.end() == prevUpMap_.find(index)) {
      if (cloud_.dataExchangeM().isFirstCouplingStep()) {
        // 如果是第一个耦合时间步，则使用初始颗粒速度
        Foam::vector Up0 = cloud_.getInitVelocity(index);
        ddtUp = (Up - Up0) / dt;
      } else {
        ddtUp = Up / dt;
      }
    } else {
      prevUp = prevUpMap_[index];
      ddtUp = (Up - prevUp) / dt;
    }
    // 计算当前时间步的 ddtUr
    ddtUr = ddtU - ddtUp;
    DDtUrHistory = Foam::vector::zero;
    if (DDtUrHistoryMap_.end() == DDtUrHistoryMap_.find(index)) {
      DDtUrHistoryMap_.insert(std::make_pair(index, std::vector<Foam::vector>()));
      DDtUrHistory = Foam::vector::zero;
    } else {
      auto& vvec = DDtUrHistoryMap_[index];
      if (1 == vvec.size()) {
        DDtUrHistory = 2 * vvec[0] / sqrt(dt);
      } else if (2 == vvec.size()) {
        DDtUrHistory = vvec[0] / sqrt(2 * dt) + vvec[1] / sqrt(dt);
      } else {
        int size = vvec.size();
        DDtUrHistory = vvec[0] / sqrt(size * dt) + vvec[size - 1] / sqrt(1 * dt);
        for (int i = 1; i < size - 1; ++i) {
          DDtUrHistory += 2 * vvec[i] / sqrt((size - i) * dt);
        }
      }
    }
    // insert ddtUr to vvec
    DDtUrHistoryMap_[index].push_back(ddtUr);
    // 计算 Basset forece
    BstForce = 0.5 * B * dt * DDtUrHistory + 2 * B * ddtUr * sqrt(dt) +
               B * initUr_ / sqrt((cloud_.dataExchangeM().couplingStep() + 1) * dt);
    if ("B" == cloud_.modelType()) {
      BstForce /= vf;
    }
    if (forceSubModel_->verbose()) {
      Pout << "Basset force = " << BstForce << endl;
    }
  }
}
#endif

}  // namespace Foam
