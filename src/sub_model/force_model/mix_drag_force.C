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

#include "./mix_drag_force.h"
#include "mpi.h"

namespace Foam {

cfdemDefineTypeName(mixDragForce);

cfdemCreateNewFunctionAdder(forceModel, mixDragForce);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
mixDragForce::mixDragForce(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      U_(cloud.globalF().U()),
      voidFraction_(cloud_.globalF().voidFraction()) {
  createForceSubModels(subPropsDict_, kUnResolved);
}

mixDragForce::~mixDragForce() {}

void mixDragForce::setForce() {
  Info << "Setting mix drag force..." << endl;
  base::MPI_Barrier();
  double dragCoefficient = 0.0;              // 阻力系数
  Foam::vector Ufluid = Foam::vector::zero;  // 背景流体速度
  Foam::vector drag = Foam::vector::zero;    // 阻力
  UInterpolator_.reset(
      interpolation<Foam::vector>::New(subPropsDict_.lookupOrDefault("UInterpolationType", word("cellPointFace")), U_)
          .ptr());
  voidFractionInterpolator_.reset(
      interpolation<Foam::scalar>::New(
          subPropsDict_.lookupOrDefault("voidfractionInterpolationType", word("cellPoint")), voidFraction_)
          .ptr());
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // 只设置 fine and middle 颗粒的阻力
    if (cloud_.checkCoarseParticle(index)) {
      continue;
    }
    // reset data
    dragCoefficient = 0.0;
    Ufluid = Foam::vector::zero;
    drag = Foam::vector::zero;
    // calculate force or drag coefficient
    setForceKernel(index, drag, Ufluid, dragCoefficient);
    // write particle data to global array
    forceSubModel_->partToArray(index, drag, Foam::vector::zero, Ufluid, dragCoefficient);
  }
  Info << "Setting mix drag force - done" << endl;
}

void mixDragForce::setForceKernel(const int index, Foam::vector& drag, Foam::vector& Ufluid, double& dragCoefficient) {
  // 颗粒中心所在网格的索引
  int findCellID = cloud_.findCellIDs()[index];
  // 获取背景流体速度
  Ufluid = getBackgroundUfluid(index, findCellID);
  // 获取背景网格空隙率
  double vf = getBackgroundVoidFraction(index, findCellID);
  // 对于每一个颗粒，只需要一个处理器计算，即颗粒中心所在的处理器
  if (findCellID >= 0) {
    double radius = cloud_.getRadius(index);               // 颗粒半径
    double diameter = 2 * radius;                          // 颗粒直径
    double nuf = forceSubModel_->nuField()[findCellID];    // 流体动力粘度
    double rho = cloud_.globalF().rhoField()[findCellID];  // 流体密度
    double pRe = 0.0;                                      // 颗粒雷诺数
    double Xi = 0.0;                                       // 模型阻力系数
    double Cd = 0.0;                                       // 流体阻力系数
    Foam::vector Up = cloud_.getVelocity(index);           // 颗粒速度
    Foam::vector Ur = Ufluid - Up;                         // 相对速度
    double magUr = mag(Ur);                                // 相对速度值
    if (magUr > 0) {
#if 0
      // 计算颗粒雷诺数
      pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
      // 计算流体阻力系数
      Cd = sqr(0.63 + 4.8 / sqrt(pRe));
      // 计算模型阻力系数
      Xi = 3.7 - 0.65 * exp(-sqr(1.5 - log10(pRe)) / 2.0);
      // 计算颗粒阻力系数
      dragCoefficient = 0.125 * Cd * rho * M_PI * diameter * diameter * pow(vf, (2 - Xi)) * magUr;
#elif 1
      // 计算颗粒雷诺数
      pRe = diameter * magUr / (nuf + Foam::SMALL);
      Cd = 24 * pow(9.06 / sqrt(pRe) + 1, 2) / (9.06 * 9.06);
      Xi = 3.7 - 0.65 * exp(-sqr(1.5 - log10(pRe)) / 2.0);
      dragCoefficient = 0.125 * Cd * rho * M_PI * diameter * diameter * pow(vf, (2 - Xi)) * magUr;
#elif 0
      // Schiller Naumann Drag
      pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
      Cd = pRe >= 1000 ? 0.44 : 24.0 * (1 + 0.15 * pow(pRe, 0.687)) / pRe;
      dragCoefficient = 0.125 * Cd * rho * M_PI * diameter * diameter * magUr;
#elif 0
      // Gidaspow drag
      pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
      if (vf > 0.8) {
        // 计算流体阻力系数
        Cd = pRe >= 1000 ? 0.44 : 24.0 * (1 + 0.15 * pow(pRe, 0.687)) / pRe;
        dragCoefficient = 0.75 * rho * vf * Cd * magUr / (diameter * Foam::pow(vf, 2.65));
      } else {
        dragCoefficient = 150 * (1 - vf) * nuf * rho / (vf * diameter * diameter) + 1.75 * magUr * rho / diameter;
      }
      dragCoefficient *= cloud_.voidFractionM().pV(radius);
#elif 0
      pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
      if (vf > 0.74) {
        double Wd = 0.0;
        if (vf <= 0.82) {
          Wd = 0.0214 / (4.0 * sqr(vf - 0.7463) + 0.0044) - 0.576;
        } else if (vf > 0.97) {
          Wd = 32.8295 * vf - 31.8295;
        } else {
          Wd = 0.0038 / (4 * sqr(vf - 0.7789) + 0.004) - 0.0101;
        }
        Cd = pRe >= 1000 ? 0.44 : (24.0 * (1 + 0.15 * pow(pRe, 0.687)) / pRe);
        dragCoefficient = 0.75 * rho * vf * Cd * magUr * Wd / diameter;
      } else {
        dragCoefficient = 150 * (1 - vf) * nuf * rho / (vf * diameter * diameter) + 1.75 * magUr * rho / diameter;
      }
      dragCoefficient *= cloud_.voidFractionM().pV(radius);
#elif 0
      pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
      double Vrs = 0.0;
      double A = 0.0;
      double B = 0.0;
      if (vf <= 0.85) {
        A = pow(vf, 4.14);
        B = 0.8 * pow(vf, 4.14);
      } else {
        A = pow(vf, 4.14);
        B = pow(vf, 2.65);
      }
      Vrs = 0.5 * (A - 0.06 * pRe + sqrt(sqr(0.06 * pRe) + 0.12 * pRe * (2 * B - A) + A * A));
      Cd = sqr(0.63 + 4.8 / sqrt(pRe / Vrs));
      dragCoefficient = 0.75 * vf * rho * Cd * magUr / (sqr(Vrs) * diameter);
      dragCoefficient *= cloud_.voidFractionM().pV(radius);
#endif
      if ("B" == cloud_.modelType()) {
        dragCoefficient /= vf;
      }
      // 计算总阻力
      drag = dragCoefficient * Ur;
    }
    if (forceSubModel_->verbose()) {
      Pout << "index = " << index << endl;
      Pout << "findCellID = " << findCellID << endl;
      Pout << "Ufluid = " << Ufluid << endl;
      Pout << "Up = " << Up << endl;
      Pout << "Ur = " << Ur << endl;
      Pout << "diameter = " << diameter << endl;
      Pout << "rho = " << rho << endl;
      Pout << "nuf = " << nuf << endl;
      Pout << "voidFraction = " << vf << endl;
      Pout << "pRe = " << pRe << endl;
      Pout << "Cd = " << Cd << endl;
      Pout << "dragCoefficient = " << dragCoefficient << endl;
      Pout << "drag (total) = " << drag << endl;
    }
  }
}

//! \brief 计算颗粒 index 处的背景流体速度
Foam::vector mixDragForce::getBackgroundUfluid(const int index, const int findCellID) const {
  // 背景流体速度
  Foam::vector Ufluid = Foam::vector::zero;
  // 对于 middle 颗粒，使用高斯核函数计算
  if (forceSubModel_->useGaussCoreFunctionRefined() && cloud_.checkMiddleParticle(index)) {
#if 1
    Ufluid = cloud_.globalF().getBackgroundUfluid(index);
#else
    if (findCellID >= 0) {
      // 获取当前颗粒的 Expanded Cell 集合
      std::unordered_set<int> set;
      double radius = cloud_.getRadius(index);
      Foam::vector particlePos = cloud_.getPosition(index);
      cloud_.voidFractionM().buildExpandedCellSet(set, findCellID, particlePos, radius, 6);
      double sumCore = 0.0;
      double core = 0.0;
      double cellV = 0.0;
      Foam::vector cellPos = Foam::vector::zero;
      for (int cellID : set) {
        if (cellID >= 0) {  // cell found
          cellPos = cloud_.mesh().C()[cellID];
          cellV = cloud_.mesh().V()[cellID];
          // 计算高斯核
          core = globalForce::GaussCore(particlePos, cellPos, radius, 6);
          // 计算累计速度
          Ufluid += voidFraction_[cellID] * U_[cellID] * core * cellV;
          // 计算累加因数
          sumCore += voidFraction_[cellID] * core * cellV;
        }
      }
      // 计算平均流体速度
      Ufluid /= sumCore;
    }
#endif
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
double mixDragForce::getBackgroundVoidFraction(const int index, const int findCellID) const {
  double vf = 1.0;
  // 对于 middle 颗粒，使用高斯核函数计算
  if (forceSubModel_->useGaussCoreFunctionRefined() && cloud_.checkMiddleParticle(index)) {
#if 1
    vf = cloud_.globalF().getBackgroundVoidFraction(index);
#else
    if (findCellID >= 0) {
      // 获取当前颗粒的 Expanded Cell 集合
      std::unordered_set<int> set;
      double radius = cloud_.getRadius(index);
      Foam::vector particlePos = cloud_.getPosition(index);
      cloud_.voidFractionM().buildExpandedCellSet(set, findCellID, particlePos, radius, 6);
      double cellV = 0.0;
      double sumPV = 0.0;
      double sumCV = 0.0;
      for (int cellID : set) {
        if (cellID >= 0) {  // cell found
          cellV = cloud_.mesh().V()[cellID];
          // 计算累加流体体积
          sumPV += voidFraction_[cellID] * cellV;
          // 计算累加网格体积
          sumCV += cellV;
        }
      }
      // 计算空隙率
      vf = sumPV / sumCV;
    }
#endif
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
