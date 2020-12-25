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

#include "./IB_void_fraction.h"

namespace Foam {

cfdemDefineTypeName(IBVoidFraction);

cfdemCreateNewFunctionAdder(voidFractionModel, IBVoidFraction);

//! \brief Constructor
IBVoidFraction::IBVoidFraction(cfdemCloud& cloud)
    : voidFractionModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      alphaMin_(subPropsDict_.lookupOrDefault<double>("alphaMin", 0.0)),
      alphaMax_(subPropsDict_.lookupOrDefault<double>("alphaMax", 1.0)) {
  // 单个颗粒覆盖最多网格数量
  maxCellsNumPerCoarseParticle_ = subPropsDict_.lookupOrDefault<int>("maxCellsNumPerCoarseParticle", 1000);
}

//! \brief Destructor
IBVoidFraction::~IBVoidFraction() {}

//! \brief 计算颗粒的体积分数场
void IBVoidFraction::setVoidFraction() {
  // reset field
  resetVolumeFraction();
  // reset particleOverMeshNumber
  base::fillTensor(cloud_.particleOverMeshNumber(), 0);
  // clear cellIDs and volumeFractions
  cloud_.pCloud().cellIDs().clear();
  cloud_.pCloud().volumeFractions().clear();
  // set volumeFraction
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // 计算index颗粒的空隙率，并设置颗粒覆盖的网格集合
    std::unique_ptr<labelHashSet> hashSetPtr(new labelHashSet);
    setVolumeFractionForSingleParticle(index, hashSetPtr);
    int meshNumber = hashSetPtr->size();
    // 检查集合中元素个数大于颗粒覆盖网格数限制
    if (meshNumber > maxCellsNumPerCoarseParticle_) {
      FatalError << __func__ << ": Big particle found " << meshNumber
                 << " cells more than permittd maximun number of cells per paticle " << maxCellsNumPerCoarseParticle_
                 << abort(FatalError);
    }
    if (meshNumber > 0) {
      // 将颗粒覆盖的当前处理器的网格数保存到 particleOverMeshNumber 中
      cloud_.particleOverMeshNumber()[index] = meshNumber;
      cloud_.pCloud().cellIDs().emplace_back(base::makeShape1(meshNumber), -1);
      cloud_.pCloud().volumeFractions().emplace_back(base::makeShape1(meshNumber), -1.0);
      for (int i = 0; i < meshNumber; ++i) {
        int cellID = hashSetPtr->toc()[i];
        // 保存颗粒覆盖的所有网格编号
        cloud_.cellIDs()[index][i] = cellID;
        // 保存 volumeFraction
        volumeFractionNext_[cellID] = volumeFractionNext_[cellID] < alphaMin_ ? alphaMin_ : volumeFractionNext_[cellID];
        volumeFractionNext_[cellID] = volumeFractionNext_[cellID] > alphaMax_ ? alphaMax_ : volumeFractionNext_[cellID];
        cloud_.volumeFractions()[index][i] = volumeFractionNext_[cellID];
      }
    } else {
      cloud_.particleOverMeshNumber()[index] = 0;
      cloud_.pCloud().cellIDs().emplace_back();
      cloud_.pCloud().volumeFractions().emplace_back();
    }
  }
}

/*!
 * \brief 设置单个颗粒的体积分数场
 * \param index 颗粒索引
 * \param hashSetPtr 哈希集合，用于保存颗粒覆盖的网格索引
 */
void IBVoidFraction::setVolumeFractionForSingleParticle(const int index,
                                                        const std::unique_ptr<labelHashSet>& hashSetPtr) {
  if (cloud_.checkPeriodicCells()) {
    FatalError << "Error: not support periodic check!" << abort(FatalError);
  }
  // 获取颗粒半径
  double radius = cloud_.getRadius(index);
  // 获取颗粒中心坐标
  Foam::vector particleCentre = cloud_.getPosition(index);
  // 获取到在当前 processor 上颗粒覆盖的某一个网格编号
  int findCellID = cloud_.findCellIDs()[index];
  if (findCellID >= 0) {  // particle centre is in domain
    // 获取网格中心坐标
    Foam::vector cellCentre = cloud_.mesh().C()[findCellID];
    // 判断网格中心是否在颗粒中
    double fc = pointInParticle(particleCentre, cellCentre, radius);
    // 计算网格的等效半径
    double corona = 0.5 * sqrt(3.0) * cbrt(cloud_.mesh().V()[findCellID]);
    // 获取网格的 corona point
    Foam::vector coronaPoint = getCoronaPointPosition(particleCentre, cellCentre, corona);
    if (pointInParticle(particleCentre, coronaPoint, radius) < 0.0) {
      // 如果 coronaPoint 在颗粒中, 则认为整个网格在颗粒中
      volumeFractionNext_[findCellID] = 0.0;
    } else {
      // 如果 coronaPoint 不在颗粒中, 则需要遍历网格的所有角点, 判断角点与网格中心是否在颗粒中
      const labelList& vertexPoints = cloud_.mesh().cellPoints()[findCellID];
      double ratio = 0.125;
      // 遍历当前网格的所有角点
      forAll(vertexPoints, i) {
        // 获取第 i 角点坐标
        vector vertexPosition = cloud_.mesh().points()[vertexPoints[i]];
        // 判断角点是否在颗粒中
        scalar fv = pointInParticle(particleCentre, vertexPosition, radius);
        if (fc < 0.0 && fv < 0.0) {
          // 网格中心在颗粒中, 角点也在颗粒中
          volumeFractionNext_[findCellID] -= ratio;
        } else if (fc < 0.0 && fv >= 0.0) {
          // 网格中心在颗粒中, 角点不在颗粒中
          // 计算角点 vertexPosition 对体积分数的影响系数 lambda
          double lambda = segmentParticleIntersection(radius, particleCentre, cellCentre, vertexPosition);
          volumeFractionNext_[findCellID] -= ratio * lambda;
        } else if (fc >= 0.0 && fv < 0.0) {
          // 网格中心不在颗粒中, 角点在颗粒中
          // 计算角点 vertexPosition 对体积分数的影响系数 lambda
          double lambda = segmentParticleIntersection(radius, particleCentre, vertexPosition, cellCentre);
          volumeFractionNext_[findCellID] -= ratio * lambda;
        }
      }  // End of loop of vertexPoints
    }
    // 颗粒中心所在网格的体积分数已经计算完成, 下面开始递归构建相邻网格
    buildLabelHashSetForVolumeFraction(findCellID, particleCentre, radius, hashSetPtr);
  }  // findCellID >= 0
}

/*!
 * \brief 构建颗粒覆盖的所有网格的哈希集合
 * \note 设置为递归函数,  通过哈希器将网格编号转换为哈希值,  并存入 set 中以便于搜索
 * \param cellID         <[in] 递归循环中要检索网格编号
 * \param particleCentre <[in] 颗粒中心位置
 * \param radius         <[in] 颗粒半径
 * \param hashSetPtr     <[in, out] 需要构建的哈希集
 */
void IBVoidFraction::buildLabelHashSetForVolumeFraction(const label cellID, const Foam::vector& particleCentre,
                                                        const double radius,
                                                        const std::unique_ptr<labelHashSet>& hashSetPtr) {
  hashSetPtr->insert(cellID);
  // 获取 cellID 网格的所有 neighbour cell 的链表
  const labelList& nc = cloud_.mesh().cellCells()[cellID];
  // 遍历 cellID 的 neighbour cell
  forAll(nc, i) {
    // 获取相邻网格索引, 以及网格中心坐标
    label neighbour = nc[i];
    if (neighbour < 0) {
      continue;
    }
    // 获取相邻网格中心坐标
    Foam::vector neighbourCentre = cloud_.mesh().C()[neighbour];
    // 判断相邻网格中心是否在颗粒中
    scalar fc = pointInParticle(particleCentre, neighbourCentre, radius);
    // 计算相邻网格的等效半径
    scalar coronaRaidus = 0.5 * sqrt(3.0) * cbrt(cloud_.mesh().V()[neighbour]);
    // 获取 corona point
    Foam::vector coronaPoint = getCoronaPointPosition(particleCentre, neighbourCentre, coronaRaidus);
    // 如果在哈希集合中没有插入 neighbour 网格
    if (false == hashSetPtr->found(neighbour)) {
      // 计算 neighbour 网格的体积分数
      if (pointInParticle(particleCentre, coronaPoint, radius) < 0.0) {
        // 如果相邻网格的 coronaPoint 在颗粒中, 则说明该网格完全被颗粒覆盖
        volumeFractionNext_[neighbour] = 0.0;
        // 以相邻网格为中心继续递归构建哈希集合
        buildLabelHashSetForVolumeFraction(neighbour, particleCentre, radius, hashSetPtr);
      } else {
        // 如果相邻网格的 coronaPoint 不在颗粒中, 则需要遍历该网格的所有角点
        // 定义单个角点对空隙率的影响率
        double ratio = 0.125;
        scalar scale = 1.0;
        // 获取 neighbour 网格的角点集合
        const labelList& vertexPoints = cloud_.mesh().cellPoints()[neighbour];
        /// 遍历网格 neighbour 的角点
        forAll(vertexPoints, j) {
          // 获取角点坐标
          Foam::vector vertexPosition = cloud_.mesh().points()[vertexPoints[j]];
          // 判断角点是否在颗粒中
          scalar fv = pointInParticle(particleCentre, vertexPosition, radius);
          if (fc < 0.0 && fv < 0.0) {  // 如果网格 neighbour 中心在颗粒中, 角点 j 也在颗粒中
            scale -= ratio;
          } else if (fc < 0.0 && fv > 0.0) {  // 如果网格 neighbour 中心在颗粒中, 角点 j 不在颗粒中
            // 计算角点对空隙率的影响系数 lambda
            scalar lambda = segmentParticleIntersection(radius, particleCentre, neighbourCentre, vertexPosition);
            scale -= lambda * ratio;
          } else if (fc > 0.0 && fv < 0.0) {  // 如果网格 neighbour 中心不在颗粒中, 角点 j 在颗粒中
            scalar lambda = segmentParticleIntersection(radius, particleCentre, vertexPosition, neighbourCentre);
            scale -= lambda * ratio;
          }
        }  // End of loop vertexPoints
        // 保证体积分数 >= 0
        scale = scale < 0.0 ? 0.0 : scale;
        scale = scale > 1.0 ? 1.0 : scale;
        if (fabs(volumeFractionNext_[neighbour] - 1.0) < Foam::SMALL) {
          // 如果 neighbour 网格的体积分数为 1.0, 则说明第一次遍历到该网格, 可以直接赋值
          volumeFractionNext_[neighbour] = scale;
        } else {
          // 如果 neighbour 网格的体积分数不为 1.0, 则说明在计算其他颗粒时候, 已经遍历到该网格
          volumeFractionNext_[neighbour] -= (1.0 - scale);
          volumeFractionNext_[neighbour] =
              volumeFractionNext_[neighbour] < alphaMin_ ? alphaMin_ : volumeFractionNext_[neighbour];
          volumeFractionNext_[neighbour] =
              volumeFractionNext_[neighbour] > alphaMax_ ? alphaMax_ : volumeFractionNext_[neighbour];
        }
        if (!(fabs(scale - 1.0) < Foam::SMALL)) {
          // 如果体积分数不为 1.0, 则说明该 neighbour 需要递归循环构建哈希集合
          buildLabelHashSetForVolumeFraction(neighbour, particleCentre, radius, hashSetPtr);
        }
      }
    }  // not found neighbour in hash set
  }    // End of loop neighbour cell
}

/*!
 * \brief 获取 Corona Point
 * \note Corona Point 是在以网格中心为中心, 半径为 corona 的球面上, 距离颗粒中心最远的点
 *       其中, 半径 corona = 0.5 * sqrt(3) * 网格等效半径
 *       网格等效半径 = pow(cellVolume, 1.0 / 3.0)
 * \param particleCentre  <[in] 指定颗粒中心
 * \param cellCentre      <[in] 指定网格中心
 * \param corona          <[in] 指定网格的等效半径
 */
Foam::vector IBVoidFraction::getCoronaPointPosition(const Foam::vector& particleCentre, const Foam::vector& cellCentre,
                                                    const scalar corona) const {
  // 计算网格中心到颗粒中心的距离
  scalar centreDist = mag(cellCentre - particleCentre);
  vector coronaPoint = cellCentre;
  if (centreDist > 0.0) {
    coronaPoint = cellCentre + (cellCentre - particleCentre) * (corona / centreDist);
    return coronaPoint;
  }
  return coronaPoint;
}

/*!
 * \brief 计算距离系数，对任意一个网格, 如果网格中心 c 在颗粒内部, 但是它的某个角点 p
 *   不在颗粒内部, 则计算 c 与 p 的连线与颗粒表面的交点 i 到网格中心 c 的距离, 即
 *   求解 x 的二元一次方程
 *   (x * (vector_p - vector_c) - vector_particle) &
 *   (x * (vector_p - vector_c) - vector_particle) == radius * radius
 *   等价于函数体中定义的: a*(x^2) - b*x + c = 0
 * \param radius         <[in] 颗粒半径
 * \param particleCentre <[in] 颗粒中心
 * \param pointInside    <[in] 网格中心
 * \param pointOutside   <[in] 网格角点
 */
double IBVoidFraction::segmentParticleIntersection(double radius, const Foam::vector& particleCentre,
                                                   const Foam::vector& pointInside,
                                                   const Foam::vector& pointOutside) const {
  // 计算方程系数 a*(x^2) - b*x + c = 0
  double a = (pointOutside - pointInside) & (pointOutside - pointInside);
  double b = 2.0 * (pointOutside - pointInside) & (pointInside - particleCentre);
  double c = ((pointInside - particleCentre) & (pointInside - particleCentre)) - radius * radius;
  double D = b * b - 4.0 * a * c;
  double lambda_1 = 0.0;
  double lambda_2 = 0.0;
  double lambda = 0.0;
  double eps = 1.0e-15;
  if (D >= 0.0) {
    // 方程有实数解, 则说明连线与球面一定有两个交点, 分别用 lambda_1 和 lambda_2 表示方程的两个实数解
    lambda_1 = (-1.0 * b + sqrt(D)) / (2.0 * a);
    lambda_2 = (-1.0 * b - sqrt(D)) / (2.0 * a);
    if (lambda_1 >= -eps && lambda_1 <= 1.0 + eps) {
      // 如果 lambda_1 在 [0, 1] 之间, 则说明这个交点在网格中心与网格角点的连线内, 则返回这个解
      lambda = lambda_1;
    } else if (lambda_2 >= -eps && lambda_2 <= 1.0 + eps) {
      // 如果 lambda_2 在 [0, 1] 之间, 则说明这个交点在网格中心与网格角点的连线内, 则返回这个解
      lambda = lambda_2;
    }
  }
  // 确保 lambda 在 [0, 1] 区间内, 而不是在 [-eps, 1 + eps] 区间内
  if (lambda > 1.0) {
    return 1;
  }
  if (lambda < 0.0) {
    return 0.0;
  }
  return lambda;
}

}  // namespace Foam
