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

#include "./mix_DiFelice_drag.h"
#include "mpi.h"

namespace Foam {

cfdemDefineTypeName(mixDiFeliceDrag);

cfdemCreateNewFunctionAdder(forceModel, mixDiFeliceDrag);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
mixDiFeliceDrag::mixDiFeliceDrag(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      velFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("velFieldName", "U").c_str()),
      voidFractionFieldName_(
          subPropsDict_.lookupOrDefault<Foam::word>("voidFractionFieldName", "voidFraction").c_str()),
      U_(cloud.mesh().lookupObject<volVectorField>(velFieldName_)),
      voidFraction_(cloud.mesh().lookupObject<volScalarField>(voidFractionFieldName_)) {
  createForceSubModels(subPropsDict_, kUnResolved);
  CHECK(false == forceSubModel_->interpolation()) << ": mixDiFeliceDrag model request interpolation == false";
}

mixDiFeliceDrag::~mixDiFeliceDrag() {}

void mixDiFeliceDrag::getMpiData(const int index, const volVectorField& field) {
  CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": ask for middle particle";
  std::unordered_set<int> set;  // 颗粒扩展网格的集合
  double radius = cloud_.getRadius(index);
  Foam::vector particlePos = cloud_.getPosition(index);
  int findExpandedCellID = cloud_.findExpandedCellIDs()[index];
  if (findExpandedCellID >= 0) {
    cloud_.voidFractionM().buildExpandedCellSet(set, findExpandedCellID, particlePos, radius, 6);
  }
}

void mixDiFeliceDrag::setForce() {
  Info << "Setting mix DiFelice drag force..." << endl;
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
    setForceKernel(index, drag, Ufluid, dragCoefficient);
    // write particle data to global array
    forceSubModel_->partToArray(index, drag, Foam::vector::zero, Ufluid, dragCoefficient);
  }
  Info << "Setting mix DiFelice drag force - done" << endl;
}

void mixDiFeliceDrag::setForceKernel(const int index, Foam::vector& drag, Foam::vector& Ufluid,
                                     double& dragCoefficient) {
  int findCellID = cloud_.findCellIDs()[index];       // 颗粒中心所在网格的索引
  Ufluid = getMpiUfluid(index, findCellID);           // 获取背景流体速度
  double vf = getMpiVoidFraction(index, findCellID);  // 获取背景网格空隙率
  base::MPI_Barrier();
  if (findCellID >= 0) {
    double radius = cloud_.getRadius(index);              // 颗粒半径
    double diameter = 2 * radius;                         // 颗粒直径
    double nuf = forceSubModel_->nuField()[findCellID];   // 流体动力粘度
    double rho = forceSubModel_->rhoField()[findCellID];  // 流体密度
    double pRe = 0.0;                                     // 颗粒雷诺数
    double Xi = 0.0;                                      // 模型阻力系数
    double Cd = 0.0;                                      // 流体阻力系数
    Foam::vector Up = cloud_.getVelocity(index);          // 颗粒速度
    Foam::vector Ur = Ufluid - Up;                        // 相对速度
    double magUr = mag(Ur);                               // 相对速度值
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
      // 这个模型在两个 test 中，结果与实验值比较吻合
      pRe = diameter * magUr * vf / (nuf + Foam::SMALL);
      Cd = 24 * pow(9.06 / sqrt(pRe) + 1, 2) / (9.06 * 9.06);
      Xi = 3.7 - 0.65 * exp(-sqr(1.5 - log10(pRe)) / 2.0);
      dragCoefficient = 0.125 * Cd * rho * M_PI * diameter * diameter * pow(vf, (2 - Xi)) * magUr;
#elif 0
      // 计算颗粒雷诺数
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
      // 计算颗粒雷诺数
      pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
      if (vf > 0.74) {
        double Wd = 0.0;
        if (vf <= 0.82) {
          Wd = 0.0214 / (4 * (vf - 0.7463) + 0.0044) - 0.576;
        } else if (vf > 0.97) {
          Wd = 32.8295 * vf - 31.8295;
        } else {
          Wd = 0.0038 / (4 * sqr(vf - 0.7789) + 0.004) - 0.0101;
        }
        Cd = pRe >= 1000 ? 0.44 : 24.0 * (1 + 0.15 * pow(pRe, 0.687)) / pRe;
        dragCoefficient = 0.75 * rho * vf * magUr * Cd * Wd / diameter;
      } else {
        dragCoefficient = 150 * (1 - vf) * nuf * rho / (vf * diameter * diameter) + 1.75 * magUr * rho / diameter;
      }
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

Foam::vector mixDiFeliceDrag::getMpiUfluid(const int index, const int findCellID) const {
  Foam::vector Ufluid = Foam::vector::zero;
  // 如果是 fine 颗粒，则需要保证 findCellID >= 0
  if (cloud_.checkFineParticle(index) && findCellID >= 0) {
    if (forceSubModel_->interpolation()) {
      // 获取颗粒中心的坐标, 将颗粒中心所在网格的空隙率和流体速度插值到颗粒中心处
      Foam::vector pos = cloud_.getPosition(index);
      Ufluid = UInterpolator_().interpolate(pos, findCellID);
    } else {
      // 不使用插值模型，直接设置颗粒中心处的值为颗粒中心所在网格的值
      Ufluid = U_[findCellID];
    }
  } else if (cloud_.checkMiddleParticle(index)) {
    Ufluid = cloud_.globalF().getBackgroundFieldValue(index, U_, voidFraction_);
    // if (findCellID >= 0) {
    //   // 获取当前颗粒的 Expanded Cell 集合
    //   std::unordered_set<int> set;
    //   double radius = cloud_.getRadius(index);
    //   Foam::vector particlePos = cloud_.getPosition(index);
    //   cloud_.voidFractionM().buildExpandedCellSet(set, findCellID, particlePos, radius, 6);
    //   double sumCore = 0.0;
    //   double core = 0.0;
    //   double cellV = 0.0;
    //   double sumPV = 0.0;
    //   double sumCV = 0.0;
    //   Foam::vector cellPos = Foam::vector::zero;
    //   for (int cellID : set) {
    //     if (cellID >= 0) {  // cell found
    //       cellPos = cloud_.mesh().C()[cellID];
    //       cellV = cloud_.mesh().V()[cellID];
    //       // 计算高斯核
    //       core = globalForce::GaussCore(particlePos, cellPos, radius, 6);
    //       // 计算累计速度
    //       Ufluid += voidFraction_[cellID] * U_[cellID] * core * cellV;
    //       // 计算累加因数
    //       sumCore += voidFraction_[cellID] * core * cellV;
    //       // 计算累加流体体积
    //       sumPV += voidFraction_[cellID] * cellV;
    //       // 计算累加网格体积
    //       sumCV += cellV;
    //     }
    //   }
    //   // 计算平均流体速度
    //   Ufluid /= sumCore;
    // }
  }
  return Ufluid;
}

double mixDiFeliceDrag::getMpiVoidFraction(const int index, const int findCellID) const {
  double vf = 1.0;
  if (cloud_.checkFineParticle(index)) {
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
  } else if (cloud_.checkMiddleParticle(index)) {
    return cloud_.globalF().getBackgroundFieldValue<false, 1, volScalarField, scalar>(index, voidFraction_,
                                                                                      voidFraction_);
    // if (findCellID >= 0) {
    //   // 获取当前颗粒的 Expanded Cell 集合
    //   std::unordered_set<int> set;
    //   double radius = cloud_.getRadius(index);
    //   Foam::vector particlePos = cloud_.getPosition(index);
    //   cloud_.voidFractionM().buildExpandedCellSet(set, findCellID, particlePos, radius, 6);
    //   double cellV = 0.0;
    //   double sumPV = 0.0;
    //   double sumCV = 0.0;
    //   for (int cellID : set) {
    //     if (cellID >= 0) {  // cell found
    //       cellV = cloud_.mesh().V()[cellID];
    //       // 计算累加流体体积
    //       sumPV += voidFraction_[cellID] * cellV;
    //       // 计算累加网格体积
    //       sumCV += cellV;
    //     }
    //   }
    //   // 计算空隙率
    //   vf = sumPV / sumCV;
    // }
  }
  return vf;
}

}  // namespace Foam
