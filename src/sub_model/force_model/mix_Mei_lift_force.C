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

#include "./mix_Mei_lift_force.h"

namespace Foam {

cfdemDefineTypeName(mixMeiLiftForce);

cfdemCreateNewFunctionAdder(forceModel, mixMeiLiftForce);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
mixMeiLiftForce::mixMeiLiftForce(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      useSecondOrderTerms_(subPropsDict_.lookupOrDefault<bool>("useSecondOrderTerms", false)),
      U_(cloud.globalF().U()),
      voidFraction_(cloud.globalF().voidFraction()),
      vorticityField_(cloud.globalF().vorticityField()) {
  createForceSubModels(subPropsDict_, kUnResolved);
}

mixMeiLiftForce::~mixMeiLiftForce() {}

void mixMeiLiftForce::setForce() {
  base::MPI_Info("Setting mix Mei lift force...", true);
  if (forceSubModel_->interpolation()) {
    // Note: not use autoPtr::reset() function for stability
    // clear interpolator before set new
    UInterpolator_.clear();
    voidFractionInterpolator_.clear();
    vorticityInterpolator_.clear();
    // set new pointer
    UInterpolator_.set(
        interpolation<Foam::vector>::New(subPropsDict_.lookupOrDefault("UInterpolationType", word("cellPoint")), U_)
            .ptr());
    voidFractionInterpolator_.set(
        interpolation<Foam::scalar>::New(
            subPropsDict_.lookupOrDefault("voidFractionInterpolationType", word("cellPoint")), voidFraction_)
            .ptr());
    vorticityInterpolator_.set(
        interpolation<Foam::vector>::New(subPropsDict_.lookupOrDefault("vorticityInterpolationType", word("cellPoint")),
                                         vorticityField_)
            .ptr());
  }
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // 只设置 fine and middle 颗粒的 lift force
    if (cloud_.checkCoarseParticle(index)) {
      continue;
    }
    // 升力
    Foam::vector mixLiftForce = Foam::vector::zero;
    setForceKernel(index, mixLiftForce);
    forceSubModel_->partToArray(index, mixLiftForce, Foam::vector::zero, Foam::vector::zero, 0);
  }
  base::MPI_Info("Setting mix Mei lift force - done", true);
}

void mixMeiLiftForce::setForceKernel(const int index, Foam::vector& mixLiftForce) {
  // 颗粒中心所在网格的索引
  int findCellID = cloud_.findCellIDs()[index];
  // 对于每一个颗粒，只需要一个处理器计算，即颗粒中心所在的处理器
  if (findCellID >= 0) {
    Foam::vector Up = cloud_.getVelocity(index);                         // 颗粒速度
    Foam::vector Ufluid = getBackgroundUfluid(index, findCellID);        // 背景流体速度
    Foam::vector vorticity = getBackgroundVorticity(index, findCellID);  // 背景流体涡量
    Foam::vector Ur = Ufluid - Up;                                       // 相对速度
    double vf = getBackgroundVoidFraction(index, findCellID);            // 背景流体空隙率
    double radius = cloud_.getRadius(index);                             // 颗粒半径
    double diameter = 2 * radius;                                        // 颗粒直径
    double rho = forceSubModel_->rhoField()[findCellID];                 // 流体密度
    double nu = forceSubModel_->nuField()[findCellID];                   // 流体动力粘度
    double magUr = mag(Ur);
    double magVorticity = mag(vorticity);
    if (magUr > 0.0 && magVorticity > 0.0) {
      double pRe = diameter * magUr / nu;                    // Loth and Dorgan (2009), Eq (3)
      double wRe = magVorticity * diameter * diameter / nu;  // Loth and Dorgan (2009), Eq (29)
      double omegaStar = magVorticity * diameter / magUr;    // Loth and Dorgan (2009), Eq (25)
      double epsilon = sqrt(magVorticity * nu) / magUr;      // McLaughlin (1991), Eq(2.7)
      double ClSaff = 0.0;                                   // Saffman lift coefficient
      double JStar = 0.0;                                    // 修正因子
      double Cl = 0.0;                                       // 修正升力系数
      // Basic model for the correction to the Saffman lift
      if (epsilon < 0.1) {  // epsilon << 1
        // McLaughlin (1991), Eq(3.27)
        // J = 32 * pi^2 * epsilon^5 * ln(1 / epsilon^2)
        // JStar = 0.443 * J
        JStar = 0.443 * 32.0 * M_PI * M_PI * pow(epsilon, 5) * log(1.0 / (epsilon * epsilon + Foam::SMALL));
      } else if (epsilon > 20.0) {  // epsilon >> 1
        // McLaughlin (1991), Eq(3.26)
        // J = 2.255 - 0.6463 / epsilon^2
        // JStar = 0.443 * J
        JStar = 0.443 * (2.255 - 0.6463 / (epsilon * epsilon));
      } else {
        // Loth and Dorgan (2009), Eq (32)
        JStar = 0.3 * (1.0 + tanh(2.5 * (log10(epsilon) + 0.191))) * (2.0 / 3.0 + tanh(6.0 * (epsilon - 0.32)));
      }
      // Saffman lift coefficient
      // Loth and Dorgan (2009), Eq (31)
      ClSaff = (12.92 / M_PI) * epsilon;
      // Loth and Dorgan (2009), Eq (32)
      Cl = JStar * ClSaff;
      // Second order terms
      if (useSecondOrderTerms_) {
        // Loth and Dorgan (2009), Eq (38)
        double omegaEq = omegaStar * (1.0 - 0.0075 * wRe) * (1 - 0.062 * sqrt(pRe) - 0.001 * pRe);
        // Loth and Dorgan (2009), Eq (34)
        double ClStar = 1.0 - (0.675 + 0.15 * (1 + tanh(0.28 * (omegaStar / 2.0 - 2.0)))) * tanh(0.18 * sqrt(pRe));
        // Loth and Dorgan (2009), Eq (39)
        Cl += omegaEq * ClStar;
      }
      // Loth and Dorgan (2009), Eq (27)
      mixLiftForce = 0.125 * M_PI * rho * Cl * magUr * (Ur ^ vorticity) * diameter * diameter / magVorticity;
      if ("B" == cloud_.modelType()) {
        mixLiftForce /= vf;
      }
      if (forceSubModel_->verbose() && 0 == index) {
        Pout << "vorticity = " << vorticity << endl;
        Pout << "wRe = " << wRe << endl;
        Pout << "epsilon = " << epsilon << endl;
        Pout << "JStar = " << JStar << endl;
        Pout << "Cl = " << Cl << endl;
        Pout << "mix lift force = " << mixLiftForce << endl;
      }
    }
  }
}

//! \brief 计算颗粒 index 处的背景流体速度
Foam::vector mixMeiLiftForce::getBackgroundUfluid(const int index, const int findCellID) const {
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

//! \brief 计算颗粒 index 处的背景流体涡量
Foam::vector mixMeiLiftForce::getBackgroundVorticity(const int index, const int findCellID) const {
  // 背景流体涡量
  Foam::vector vorticity = Foam::vector::zero;
  // 对于 middle 颗粒，使用高斯核函数计算
  if (forceSubModel_->useGaussCoreFunctionRefined() && cloud_.checkMiddleParticle(index)) {
    vorticity = cloud_.globalF().getBackgroundVorticity(index);
  } else if (!cloud_.checkCoarseParticle(index) && findCellID >= 0) {
    if (forceSubModel_->interpolation()) {
      // 获取颗粒中心的坐标, 将颗粒中心所在网格的流体涡量插值到颗粒中心处
      Foam::vector pos = cloud_.getPosition(index);
      vorticity = vorticityInterpolator_().interpolate(pos, findCellID);
    } else {
      // 不使用插值模型，直接设置颗粒中心处的值为颗粒中心所在网格的值
      vorticity = vorticityField_[findCellID];
    }
  }
  return vorticity;
}

//! \brief 计算颗粒 index 处的背景流体空隙率
double mixMeiLiftForce::getBackgroundVoidFraction(const int index, const int findCellID) const {
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
