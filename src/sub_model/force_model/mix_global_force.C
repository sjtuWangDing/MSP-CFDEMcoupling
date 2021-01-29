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

#include "./mix_Basset_force.h"
#include "./mix_Mei_lift_force.h"
#include "./mix_global_force.h"
#include "./mix_virtual_mass_force.h"
#include "cfdem_tools/cfdem_tools.h"
#include "sub_model/void_fraction_model/void_fraction_model.h"

namespace Foam {

cfdemDefineTypeName(mixGlobalForce);

cfdemCreateNewFunctionAdder(globalForce, mixGlobalForce);

//! \brief Constructor
mixGlobalForce::mixGlobalForce(cfdemCloud& cloud) : globalForce(cloud) {}

//! \brief Destructor
mixGlobalForce::~mixGlobalForce() {}

//! \brief 每一次耦合中，在 set force 前执行
void mixGlobalForce::initBeforeSetForce() {
  // init data - Important !!!
  // nerver forget !!!
  expandedCellMap_.clear();
  backgroundUfluidMap_.clear();
  backgroundVoidFractionMap_.clear();
  backgroundDDtUMap_.clear();
  backgroundVorticityMap_.clear();

  bool isUsedVirtualMassForce = cfdemTools::isUsedForceModel(cloud_, mixVirtualMassForce::cTypeName());
  bool isUsedBassetForce = cfdemTools::isUsedForceModel(cloud_, mixBassetForce::cTypeName());
  // 如果使用了 BassetForce or virtualMassForce，则需要计算 ddtU
  if (isUsedVirtualMassForce || isUsedBassetForce) {
    // 计算 ddtU field
    ddtU_ = fvc::ddt(U_) + fvc::div(phi_, U_);
  }
  // 如果使用了升力，则需要计算 vorticityField
  bool isUsedMixMeiLiftForce = cfdemTools::isUsedForceModel(cloud_, mixMeiLiftForce::cTypeName());
  if (isUsedMixMeiLiftForce) {
    // 计算 vorticityField
    vorticityField_ = fvc::curl(U_);
  }
  // reset data
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkMiddleParticle(index)) {
      double radius = cloud_.getRadius(index);                       // 颗粒半径
      Foam::vector particlePos = cloud_.getPosition(index);          // 颗粒中心坐标
      int findExpandedCellID = cloud_.findExpandedCellIDs()[index];  // 扩展网格ID
      // 计算颗粒覆盖的扩展网格集合
      expandedCellMap_.insert(std::make_pair(index, std::unordered_set<int>()));
      if (findExpandedCellID >= 0) {
        cloud_.voidFractionM().buildExpandedCellSet(expandedCellMap_[index], findExpandedCellID, particlePos, radius,
                                                    cloud_.expandedCellScale());
      }
      // 计算背景流体速度
      backgroundUfluidMap_.insert(std::make_pair(index, getBackgroundFieldValue(index, U_)));
      // 计算背景流体空隙率
      backgroundVoidFractionMap_.insert(
          std::make_pair(index, getBackgroundFieldValue<false, 1, volScalarField, scalar>(index, voidFraction_)));
      // 计算背景流体的 ddtU
      if (isUsedVirtualMassForce || isUsedBassetForce) {
        backgroundDDtUMap_.insert(std::make_pair(index, getBackgroundFieldValue(index, ddtU_)));
      }
      // 计算背景流体的涡量
      if (isUsedMixMeiLiftForce) {
        backgroundVorticityMap_.insert(std::make_pair(index, getBackgroundFieldValue(index, vorticityField_)));
      }
    }
  }
  base::MPI_Info("mixGlobalForce: initBeforeSetForce - done", verbose_);
}

//! \brief 每一次耦合中，在 set force 后执行
void mixGlobalForce::endAfterSetForce() {
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // (1) 更新颗粒速度，所有处理器都必须更新，当颗粒从一个计算域运动到另一个计算域，新的计算域必须有颗粒的 prevU
    updatePrevParticleVelMap(index);
    // (2) 更新 ddtUrHistoryMap_，所有处理器都必须更新
    if (cloud_.checkMiddleParticle(index) && cfdemTools::isUsedForceModel(cloud_, mixBassetForce::cTypeName())) {
      updateDDtUrHistoryMap(index);
    }
  }
  base::MPI_Info("mixGlobalForce: endAfterSetForce - done", verbose_);
}

//! \brief 获取颗粒处背景流体速度
Foam::vector mixGlobalForce::getBackgroundUfluid(const int index) const {
  CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
  auto iter = backgroundUfluidMap_.find(index);
  if (backgroundUfluidMap_.end() != iter) {
    return iter->second;
  }
  return Foam::vector::zero;
}

//! \brief 获取颗粒处背景空隙率
double mixGlobalForce::getBackgroundVoidFraction(const int index) const {
  CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
  auto iter = backgroundVoidFractionMap_.find(index);
  if (backgroundVoidFractionMap_.end() != iter) {
    return iter->second;
  }
  return 1.0;
}

//! \brief 获取颗粒处背景流体的 ddtU
Foam::vector mixGlobalForce::getBackgroundDDtU(const int index) const {
  CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
  auto iter = backgroundDDtUMap_.find(index);
  if (backgroundDDtUMap_.end() != iter) {
    return iter->second;
  }
  return Foam::vector::zero;
}

//! \brief 获取颗粒处背景流体的涡量
Foam::vector mixGlobalForce::getBackgroundVorticity(const int index) const {
  CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
  auto iter = backgroundVorticityMap_.find(index);
  if (backgroundVorticityMap_.end() != iter) {
    return iter->second;
  }
  return Foam::vector::zero;
}

}  // namespace Foam