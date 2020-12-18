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
  demoModel
\*---------------------------------------------------------------------------*/

#include "./void_fraction_model.h"

namespace Foam {

cfdemDefineTypeName(voidFractionModel);

cfdemDefineNewFunctionMap(voidFractionModel);

cfdemDefineConstructNewFunctionMap(voidFractionModel);

cfdemDefineDestroyNewFunctionMap(voidFractionModel);

cfdmeDefineBaseTypeNew(autoPtr, voidFractionModel, (cfdemCloud & cloud, const dictionary& dict), dict, (cloud));

//! \brief Constructor
voidFractionModel::voidFractionModel(cfdemCloud& cloud)
    : cloud_(cloud),
      voidFractionPrev_(IOobject("voidFractionPrev", cloud.mesh().time().timeName(), cloud.mesh(),
                                 IOobject::READ_IF_PRESENT,  // or MUST_READ,
                                 IOobject::AUTO_WRITE),
                        // cloud.mesh().lookupObject<volScalarField>("voidFraction")
                        cloud.mesh(), dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)),
      voidFractionNext_(IOobject("voidFractionNext", cloud.mesh().time().timeName(), cloud.mesh(),
                                 IOobject::READ_IF_PRESENT,  // or MUST_READ,
                                 IOobject::AUTO_WRITE),
                        // cloud.mesh().lookupObject<volScalarField>("voidFraction")
                        cloud.mesh(), dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)),
      volumeFractionPrev_(IOobject("volumeFractionPrev", cloud.mesh().time().timeName(), cloud.mesh(),
                                   IOobject::READ_IF_PRESENT,  // or MUST_READ,
                                   IOobject::AUTO_WRITE),
                          // cloud.mesh().lookupObject<volScalarField>("volumeFraction")
                          cloud.mesh(), dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)),
      volumeFractionNext_(IOobject("volumeFractionNext", cloud.mesh().time().timeName(), cloud.mesh(),
                                   IOobject::READ_IF_PRESENT,  // or MUST_READ,
                                   IOobject::AUTO_WRITE),
                          // cloud.mesh().lookupObject<volScalarField>("volumeFraction")
                          cloud.mesh(), dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)),
      maxCellsNumPerFineParticle_(1),
      maxCellsNumPerMiddleParticle_(1),
      maxCellsNumPerCoarseParticle_(1) {}

//! \brief Destructor
voidFractionModel::~voidFractionModel() {}

/*!
 * \brief 构建颗粒覆盖的所有网格的哈希集合
 * \note 设置为递归函数，通过哈希器将网格编号转换为哈希值，并存入 set 中以便于搜索
 * \param hashSett    <[in, out] 需要构建的哈希集
 * \param cellID      <[in] 递归循环中要检索网格编号
 * \param particlePos <[in] 颗粒中心位置
 * \param radius      <[in] 颗粒半径
 * \param scale       <[in] 颗粒半径扩大系数
 */
void voidFractionModel::buildLabelHashSetForCoveredCell(labelHashSet& hashSett, const int cellID,
                                                        const Foam::vector& particlePos, const double radius,
                                                        const double scale) const {
  // 如果搜索到网格的边界，则结束搜索
  if (cellID < 0) {
    return;
  }
  // 将 cellID 插入到哈希集合中
  hashSett.insert(cellID);
  // 获取网格中心坐标
  Foam::vector neighbourPos = Foam::vector::zero;
  // 获取 cellID 网格的所有 neighbour cell 的链表
  const labelList& nc = cloud_.mesh().cellCells()[cellID];
  // 遍历链表
  for (int i = 0; i < nc.size(); ++i) {
    int neighbour = nc[i];
    neighbourPos = cloud_.mesh().C()[neighbour];
    if (false == hashSett.found(neighbour) && pointInParticle(neighbourPos, particlePos, radius, scale) < 0.0) {
      // 如果 neighbour 网格中心在颗粒中, 并且在哈希集合中没有 neighbour 网格的索引
      // 以 neighbour 为中心递归构建哈希集合
      buildLabelHashSetForCoveredCell(hashSett, neighbour, particlePos, radius, scale);
    }
  }  // end loop of neighbour list
}

//! \brief 计算颗粒尺寸与其周围网格平均尺寸的比值, 并将颗粒索引按照颗粒尺寸归类
void voidFractionModel::getDimensionRatios(const base::CITensor1& findCellIDs, const base::CDTensor1& dimensionRatios,
                                           const base::CITensor1& particleOverMeshNumber,
                                           const double scale /* = 1.0 */) const {
  std::fill_n(dimensionRatios.ptr(), dimensionRatios.mSize(), 0.0);
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // init
    particleOverMeshNumber[index] = 0;
    dimensionRatios[index] = 0.0;
    int findCellID = findCellIDs[index];
    if (findCellID < 0) {
      continue;
    }
    // 获取颗粒半径
    double radius = cloud_.getRadius(index);
    // 获取颗粒中心位置
    Foam::vector particlePos = cloud_.getPosition(index);
    // 定义哈希集合
    labelHashSet initHashSett;
    // 构建初始化哈希集合
    buildLabelHashSetForCoveredCell(initHashSett, findCellID, particlePos, radius, scale);
    // 计算颗粒周围网格的平均尺寸
    // Pout << "particleCenterCellID: " << findCellID << ", " << initHashSett.size() << endl;
    particleOverMeshNumber[index] = initHashSett.size();
  }  // end loop of particles
}

}  // namespace Foam
