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

#include "./mix_grad_p_force.h"

namespace Foam {

cfdemDefineTypeName(mixGradPForce);

cfdemCreateNewFunctionAdder(forceModel, mixGradPForce);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
mixGradPForce::mixGradPForce(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      p_(cloud.globalF().p()),
      gradPField_(cloud_.globalF().gradPField()) {
  createForceSubModels(subPropsDict_, kUnResolved);
  if (cloud_.modelType() == "B") {
    FatalError << "using mixGradPForce with model type B is not valid\n" << abort(FatalError);
  }
  if (!forceSubModel_->treatForceExplicitInMomEquation()) {
    FatalError << "mixGradPForce model need treatForceExplicitInMomEquation = true\n" << abort(FatalError);
  }
  if (cloud_.modelType() == "A" && forceSubModel_->treatForceBothCFDAndDEM()) {
    FatalError << "mixGradPForce with model type A requires treatForceBothCFDAndDEM = false\n" << abort(FatalError);
  }
  if (cloud_.modelType() == "Bfull" && !forceSubModel_->treatForceBothCFDAndDEM()) {
    FatalError << "mixGradPForce with model type Bfull requires treatForceBothCFDAndDEM = true\n" << abort(FatalError);
  }
}

mixGradPForce::~mixGradPForce() {}

void mixGradPForce::setForce() {
  base::MPI_Info("Setting mix grad p force...", true);
  if (forceSubModel_->interpolation()) {
    gradPInterpolator_.clear();
    gradPInterpolator_.set(interpolation<Foam::vector>::New(
                               subPropsDict_.lookupOrDefault("gradPInterpolationType", word("cellPoint")), gradPField_)
                               .ptr());
  }
  std::once_flag onceOp;
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // 只设置 fine and middle 颗粒
    if (cloud_.checkCoarseParticle(index)) {
      continue;
    }
    // 压力梯度力
    Foam::vector mixGradPF = Foam::vector::zero;
    setForceKernel(index, mixGradPF, onceOp);
    // write particle data to global array
    forceSubModel_->partToArray(index, mixGradPF, Foam::vector::zero, Foam::vector::zero, 0.0);
  }
  base::MPI_Info("Setting mix grad p force - done", true);
}

void mixGradPForce::setForceKernel(const int index, Foam::vector& mixGradPF, std::once_flag& onceOp) {
  // 颗粒中心所在网格的索引
  int findCellID = cloud_.findCellIDs()[index];
  // 对于每一个颗粒，只需要一个处理器计算，即颗粒中心所在的处理器
  if (findCellID >= 0) {
    double rho = forceSubModel_->rhoField()[findCellID];  // 流体密度
    double radius = cloud_.getRadius(index);              // 颗粒半径
    Foam::vector gradP = getBackgroundGradP(index, findCellID);
    mixGradPF = -gradP * rho * cloud_.voidFractionM().pV(radius);
    std::call_once(onceOp, [&]() {
      if (forceSubModel_->verbose() && mag(mixGradPF) > Foam::SMALL) {
        Pout << "mix grad p = " << gradP << endl;
        Pout << "mix grad p force = " << mixGradPF << endl;
      }
    });
  }
}

Foam::vector mixGradPForce::getBackgroundGradP(const int index, const int findCellID) const {
  // 背景流体 gradP
  Foam::vector gradP = Foam::vector::zero;
  // 对于 middle 颗粒，使用高斯核函数计算
  if (forceSubModel_->useGaussCoreFunctionRefined() && cloud_.checkMiddleParticle(index)) {
    gradP = cloud_.globalF().getBackgroundGradP(index);
  } else if (!cloud_.checkCoarseParticle(index) && findCellID >= 0) {
    if (forceSubModel_->interpolation()) {
      Foam::vector pos = cloud_.getPosition(index);
      gradP = gradPInterpolator_().interpolate(pos, findCellID);
    } else {
      // 不使用插值模型，直接设置颗粒中心处的值为颗粒中心所在网格的值
      gradP = gradPField_[findCellID];
    }
  }
  return gradP;
}

}  // namespace Foam
