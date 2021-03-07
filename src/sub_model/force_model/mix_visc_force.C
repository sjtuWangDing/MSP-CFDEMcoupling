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

#include "./mix_visc_force.h"

namespace Foam {

cfdemDefineTypeName(mixViscForce);

cfdemCreateNewFunctionAdder(forceModel, mixViscForce);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
mixViscForce::mixViscForce(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      divTauField_(cloud_.globalF().divTauField()) {
  createForceSubModels(subPropsDict_, kUnResolved);
  if (cloud_.modelType() == "B") {
    FatalError << "using mixViscForce with model type B is not valid\n" << abort(FatalError);
  }
  // if (!forceSubModel_->treatForceExplicitInMomEquation()) {
  //   FatalError << "mixViscForce model need treatForceExplicitInMomEquation = true\n" << abort(FatalError);
  // }
  if (cloud_.modelType() == "A" && forceSubModel_->treatForceBothCFDAndDEM()) {
    FatalError << "mixViscForce with model type A requires treatForceBothCFDAndDEM = false\n" << abort(FatalError);
  }
  if (cloud_.modelType() == "Bfull" && !forceSubModel_->treatForceBothCFDAndDEM()) {
    FatalError << "mixViscForce with model type Bfull requires treatForceBothCFDAndDEM = true\n" << abort(FatalError);
  }
}

mixViscForce::~mixViscForce() {}

void mixViscForce::setForce() {
  base::MPI_Info("Setting mix visc force...", true);
  if (forceSubModel_->interpolation()) {
    divTauInterpolator_.clear();
    divTauInterpolator_.set(
        interpolation<Foam::vector>::New(subPropsDict_.lookupOrDefault("divTauInterpolationType", word("cellPoint")),
                                         divTauField_)
            .ptr());
  }
  std::once_flag onceOp;
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // 只设置 fine and middle 颗粒
    if (cloud_.checkCoarseParticle(index)) {
      continue;
    }
    // 粘性应力
    Foam::vector viscForce = Foam::vector::zero;
    setForceKernel(index, viscForce, onceOp);
    forceSubModel_->partToArray(index, viscForce, Foam::vector::zero, Foam::vector::zero, 0.0);
  }
  base::MPI_Info("Setting mix visc force - done", true);
}

void mixViscForce::setForceKernel(const int index, Foam::vector& viscForce, std::once_flag& onceOp) {
  // 颗粒中心所在网格的索引
  int findCellID = cloud_.findCellIDs()[index];
  // 对于每一个颗粒，只需要一个处理器计算，即颗粒中心所在的处理器
  if (findCellID >= 0) {
    double radius = cloud_.getRadius(index);  // 颗粒半径
    Foam::vector divTau = getBackgroundDivTau(index, findCellID);
    viscForce = -divTau * cloud_.voidFractionM().pV(radius);
    std::call_once(onceOp, [&]() {
      if (forceSubModel_->verbose() && mag(viscForce) > Foam::SMALL) {
        Pout << "mix divTau = " << divTau << endl;
        Pout << "mix visc force = " << viscForce << endl;
      }
    });
  }
}

Foam::vector mixViscForce::getBackgroundDivTau(const int index, const int findCellID) const {
  // 背景流体 divTau
  Foam::vector divTau = Foam::vector::zero;
  // 对于 middle 颗粒，使用高斯核函数计算
  if (forceSubModel_->useGaussCoreFunctionRefined() && cloud_.checkMiddleParticle(index)) {
    divTau = cloud_.globalF().getBackgroundDivTau(index);
  } else if (!cloud_.checkCoarseParticle(index) && findCellID >= 0) {
    if (forceSubModel_->interpolation()) {
      Foam::vector pos = cloud_.getPosition(index);
      divTau = divTauInterpolator_().interpolate(pos, findCellID);
    } else {
      // 不使用插值模型，直接设置颗粒中心处的值为颗粒中心所在网格的值
      divTau = divTauField_[findCellID];
    }
  }
  return divTau;
}

}  // namespace Foam
