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
#include "./divided_void_fraction.h"
#include "./mix_void_fraction.h"
#include "sub_model/locate_model/locate_model.h"

namespace Foam {

cfdemDefineTypeName(mixVoidFraction);

cfdemCreateNewFunctionAdder(voidFractionModel, mixVoidFraction);

//! \brief Constructor
mixVoidFraction::mixVoidFraction(cfdemCloud& cloud)
    : voidFractionModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      verbose_(false),
      alphaMin_(0.0),
      tooMuch_(0.0) {
  weight_ = subPropsDict_.lookupOrDefault<double>("weight", 1.0);
  porosity_ = subPropsDict_.lookupOrDefault<double>("porosity", 1.0);
  verbose_ = subPropsDict_.lookupOrDefault<bool>("verbose", false);
  alphaMin_ = subPropsDict_.lookupOrDefault<double>("alphaMin", 0.0);
  // 单个颗粒覆盖最多网格数量
  maxCellsNumPerCoarseParticle_ = subPropsDict_.lookupOrDefault<int>("maxCellsNumPerCoarseParticle", 1000);
}

//! \brief Destructor
mixVoidFraction::~mixVoidFraction() {}

//! \brief 输出空隙率相关信息
void mixVoidFraction::printVoidFractionInfo() const {
  int nProcs = 0, id = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  base::MPI_Barrier();
  for (int i = 0; i < nProcs; ++i) {
    if (id == i) {
      Pout << typeName().c_str() << ":" << endl;
      for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
        Pout << "  findCellIDs[" << index << "]: " << cloud_.findCellIDs()[index] << endl;
        Pout << "  particleOverMeshNumber[" << index << "]: " << cloud_.particleOverMeshNumber()[index] << endl;
        // Pout << "  cellIDs[" << index << "]: ";
        // for (int j = 0; j < cloud_.particleOverMeshNumber()[index]; ++j) {
        //   Pout << cloud_.cellIDs()[index][j] << ", ";
        // }
        // Pout << endl;
        if (cloud_.checkFineParticle(index) || cloud_.checkMiddleParticle(index)) {
          Pout << "  voidFractions[" << index << "]: ";
          for (int j = 0; j < cloud_.particleOverMeshNumber()[index]; ++j) {
            Pout << cloud_.voidFractions()[index][j] << ", ";
          }
          Pout << endl;
          Pout << "  particleWeights[" << index << "]: ";
          for (int j = 0; j < cloud_.particleOverMeshNumber()[index]; ++j) {
            Pout << cloud_.particleWeights()[index][j] << ", ";
          }
          Pout << endl;
          Pout << "  particleVolumes[" << index << "]: ";
          for (int j = 0; j < cloud_.particleOverMeshNumber()[index]; ++j) {
            Pout << cloud_.particleVolumes()[index][j] << ", ";
          }
          Pout << endl;
        }
      }
    }
    base::MPI_Barrier();
  }  // End of procs loop
}

void mixVoidFraction::setVoidFraction() {
  // reset field
  resetVoidFraction();
  resetVolumeFraction();
  // reset particleOverMeshNumber
  base::fillTensor(cloud_.particleOverMeshNumber(), 0);
  // clear before set
  cloud_.pCloud().cellIDs().clear();
  cloud_.pCloud().voidFractions().clear();
  cloud_.pCloud().volumeFractions().clear();
  cloud_.pCloud().particleWeights().clear();
  cloud_.pCloud().particleVolumes().clear();
  // define data buffer
  // parMaps for fine and middle particles and parSets for coarse particles
  std::unordered_map<int, std::unordered_map<int, Foam::vector>> parMaps;
  std::unordered_map<int, std::unique_ptr<labelHashSet>> parSets;
  // set voidFraction and volumeFraction
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkCoarseParticle(index)) {
      // coarse particle
      std::unique_ptr<labelHashSet> hashSetPtr(new labelHashSet);
      setVolumeFractionForSingleParticle(index, hashSetPtr);
      parSets.insert(std::make_pair(index, std::move(hashSetPtr)));
    } else {
      // find and middle particle
      int findCellID = cloud_.findCellIDs()[index];
      parMaps.insert(std::make_pair(index, std::unordered_map<int, Foam::vector>()));
      if (findCellID >= 0) {
        setVoidFractionForSingleParticle(index, findCellID, parMaps[index]);
      }
    }
    // no matter what kind of particle, need to construct default tensor
    cloud_.pCloud().cellIDs().emplace_back();
    cloud_.pCloud().voidFractions().emplace_back();
    cloud_.pCloud().volumeFractions().emplace_back();
    cloud_.pCloud().particleWeights().emplace_back();
    cloud_.pCloud().particleVolumes().emplace_back();
  }  // End loop of all particles
  voidFractionNext_.correctBoundaryConditions();
  volumeFractionNext_.correctBoundaryConditions();
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkCoarseParticle(index)) {
      const auto& hashSetPtr = parSets[index];
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
        // allocate memory for coarse particle
        cloud_.pCloud().cellIDs()[index] = std::move(base::CITensor1(base::makeShape1(meshNumber), -1));
        cloud_.pCloud().volumeFractions()[index] = std::move(base::CDTensor1(base::makeShape1(meshNumber), 0.0));
        for (int i = 0; i < meshNumber; ++i) {
          int cellID = hashSetPtr->toc()[i];
          // 保存颗粒覆盖的所有网格编号
          cloud_.cellIDs()[index][i] = cellID;
          // 保存 volumeFraction
          cloud_.volumeFractions()[index][i] = volumeFractionNext_[cellID];
        }
      }
    } else {
      const auto& umap = parMaps[index];
      int meshNumber = umap.size();
      if (meshNumber > 0) {
        // realloc memory according to meshNumber
        cloud_.particleOverMeshNumber()[index] = meshNumber;
        cloud_.pCloud().cellIDs()[index] = std::move(base::CITensor1(base::makeShape1(meshNumber), -1));
        cloud_.pCloud().voidFractions()[index] = std::move(base::CDTensor1(base::makeShape1(meshNumber), 0.0));
        cloud_.pCloud().particleWeights()[index] = std::move(base::CDTensor1(base::makeShape1(meshNumber), 0.0));
        cloud_.pCloud().particleVolumes()[index] = std::move(base::CDTensor1(base::makeShape1(meshNumber), 0.0));
        // 遍历 umap
        auto iter = umap.begin();
        for (int i = 0; i < meshNumber; ++i, ++iter) {
          int subCellID = iter->first;
          const Foam::vector& data = iter->second;
          // 保存颗粒覆盖的所有网格编号
          cloud_.cellIDs()[index][i] = subCellID;
          // 保存 voidFractions
          cloud_.voidFractions()[index][i] = voidFractionNext_[subCellID];
          // 保存 particleWeights
          cloud_.particleWeights()[index][i] = data[0];
          // 保存 particleVolumes
          cloud_.particleVolumes()[index][i] = data[1];
        }
      }
    }
  }  // End loop of all particles
  base::MPI_Info("mixVoidFraction: setVoidFraction - done", verbose_);
}

/*!
 * \brief 设置单个颗粒的体积分数场
 * \param index 颗粒索引
 * \param hashSetPtr 哈希集合，用于保存颗粒覆盖的网格索引
 */
void mixVoidFraction::setVolumeFractionForSingleParticle(const int index,
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
    double fc = pointInParticle(cellCentre, particleCentre, radius);
    // 计算网格的等效半径
    double corona = 0.5 * sqrt(3.0) * cbrt(cloud_.mesh().V()[findCellID]);
    // 获取网格的 corona point
    Foam::vector coronaPoint = IBVoidFraction::getCoronaPointPosition(particleCentre, cellCentre, corona);
    if (pointInParticle(coronaPoint, particleCentre, radius) < 0.0) {
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
        scalar fv = pointInParticle(vertexPosition, particleCentre, radius);
        if (fc < 0.0 && fv < 0.0) {
          // 网格中心在颗粒中, 角点也在颗粒中
          volumeFractionNext_[findCellID] -= ratio;
        } else if (fc < 0.0 && fv >= 0.0) {
          // 网格中心在颗粒中, 角点不在颗粒中
          // 计算角点 vertexPosition 对体积分数的影响系数 lambda
          double lambda =
              IBVoidFraction::segmentParticleIntersection(radius, particleCentre, cellCentre, vertexPosition);
          volumeFractionNext_[findCellID] -= ratio * lambda;
        } else if (fc >= 0.0 && fv < 0.0) {
          // 网格中心不在颗粒中, 角点在颗粒中
          // 计算角点 vertexPosition 对体积分数的影响系数 lambda
          double lambda =
              IBVoidFraction::segmentParticleIntersection(radius, particleCentre, vertexPosition, cellCentre);
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
 * \note 设置为递归函数,  通过哈希器将网格编号转换为哈希值, 并存入 set 中以便于搜索
 * \param cellID         <[in] 递归循环中要检索网格编号
 * \param particleCentre <[in] 颗粒中心位置
 * \param radius         <[in] 颗粒半径
 * \param hashSetPtr     <[in, out] 需要构建的哈希集
 */
void mixVoidFraction::buildLabelHashSetForVolumeFraction(const label cellID, const Foam::vector& particleCentre,
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
    scalar fc = pointInParticle(neighbourCentre, particleCentre, radius);
    // 计算相邻网格的等效半径
    scalar coronaRaidus = 0.5 * sqrt(3.0) * cbrt(cloud_.mesh().V()[neighbour]);
    // 获取 corona point
    Foam::vector coronaPoint = IBVoidFraction::getCoronaPointPosition(particleCentre, neighbourCentre, coronaRaidus);
    // 如果在哈希集合中没有插入 neighbour 网格
    if (false == hashSetPtr->found(neighbour)) {
      // 计算 neighbour 网格的体积分数
      if (pointInParticle(coronaPoint, particleCentre, radius) < 0.0) {
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
          scalar fv = pointInParticle(vertexPosition, particleCentre, radius);
          if (fc < 0.0 && fv < 0.0) {  // 如果网格 neighbour 中心在颗粒中, 角点 j 也在颗粒中
            scale -= ratio;
          } else if (fc < 0.0 && fv > 0.0) {  // 如果网格 neighbour 中心在颗粒中, 角点 j 不在颗粒中
            // 计算角点对空隙率的影响系数 lambda
            scalar lambda =
                IBVoidFraction::segmentParticleIntersection(radius, particleCentre, neighbourCentre, vertexPosition);
            scale -= lambda * ratio;
          } else if (fc > 0.0 && fv < 0.0) {  // 如果网格 neighbour 中心不在颗粒中, 角点 j 在颗粒中
            scalar lambda =
                IBVoidFraction::segmentParticleIntersection(radius, particleCentre, vertexPosition, neighbourCentre);
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
        }
        if (!(fabs(scale - 1.0) < Foam::SMALL)) {
          // 如果体积分数不为 1.0, 则说明该 neighbour 需要递归循环构建哈希集合
          buildLabelHashSetForVolumeFraction(neighbour, particleCentre, radius, hashSetPtr);
        }
      }
    }  // not found neighbour in hash set
  }    // End of loop neighbour cell
}

//! \brief 设置索引为 index 的单个颗粒的空隙率
//! \note used for fine particle
void mixVoidFraction::setVoidFractionForSingleParticle(const int index, const int findCellID,
                                                       std::unordered_map<int, Foam::vector>& parMap) {
  parMap.clear();
  int numberOfMarkerPoints = dividedVoidFraction::numberOfMarkerPoints();             // 标志点个数
  int findSubCellNum = 0;                                                             // 有效标记点的个数
  double radius = cloud_.getRadius(index);                                            // 颗粒半径
  double scaleVol = weight();                                                         // 颗粒体积因子
  double scaleRadius = cbrt(porosity());                                              // 颗粒半径因子
  double avgPVol = pV(radius, scaleVol) / static_cast<double>(numberOfMarkerPoints);  // 获取平均颗粒体积（1 / 29）
  Foam::vector particleCentre = cloud_.getPosition(index);                            // 颗粒中心坐标
  Foam::vector subPosition = Foam::vector::zero;                                      // 标志点的绝对位置
  // 先计算体积, 然后再用 scaleRadius 乘以半径, 以保证在计算空隙率时候颗粒尺寸不变
  radius *= scaleRadius;
  // 遍历所有标志点，除了颗粒中心处的标志点
  for (int i = 1; i < numberOfMarkerPoints; ++i) {
    // 获取标志点的绝对位置
    subPosition = particleCentre + radius * dividedVoidFraction::offsets()[i];
    if (cloud_.checkPeriodicCells()) {
      FatalError << "Error: not support periodic check!" << abort(FatalError);
    }
    // 根据修正后的标志点坐标定位标志点所在网格的索引
    int subCellId = cloud_.locateM().findSingleCell(subPosition, findCellID);
    // 如果标志点在 domain 中, 则更新空隙率
    if (subCellId >= 0) {
      findSubCellNum += 1;
      // 获取 subCellId 网格体积
      double subCellVol = cloud_.mesh().V()[subCellId];
      double newAlpha = voidFractionNext_[subCellId] - avgPVol / subCellVol;
      // 更新空隙率场
      if (newAlpha > alphaMin_) {
        voidFractionNext_[subCellId] = newAlpha;
      } else {
        // 如果空隙率低于最低空隙率, 则直接赋值为 alphaMin_, 并累加损失的颗粒体积
        voidFractionNext_[subCellId] = alphaMin_;
        tooMuch_ += (alphaMin_ - newAlpha) * subCellVol;
      }
      auto iter = parMap.find(subCellId);
      if (parMap.end() == iter) {
        // 需要创建新的索引
        cloud_.particleOverMeshNumber()[index] += 1;
        parMap.insert(std::make_pair(subCellId, Foam::vector::zero));
      }
      parMap[subCellId][0] += 1.0 / static_cast<double>(numberOfMarkerPoints);  // add particleWeights
      parMap[subCellId][1] += avgPVol;                                          // add particleVolumes
      parMap[subCellId][2] += avgPVol;                                          // add particleVs
    }
  }  // End of marker point loop
  // 处理颗粒中心处的标志点, 将丢失的标志点对空隙率的影响作用在颗粒中心处的标志点上
  // set source for particle center
  double findCellVol = cloud_.mesh().V()[findCellID];
  double centreWeight = (numberOfMarkerPoints - findSubCellNum) * (1.0 / numberOfMarkerPoints);
  double newAlpha = voidFractionNext_[findCellID] - pV(radius, scaleVol) * centreWeight / findCellVol;
  // update voidFraction for findCellID
  if (newAlpha > alphaMin_) {
    voidFractionNext_[findCellID] = newAlpha;
  } else {
    voidFractionNext_[findCellID] = alphaMin_;
    tooMuch_ += (alphaMin_ - newAlpha) * findCellVol;
  }
  auto iter = parMap.find(findCellID);
  if (parMap.end() == iter) {
    // 需要创建新的索引
    cloud_.particleOverMeshNumber()[index] += 1;
    parMap.insert(std::make_pair(findCellID, Foam::vector::zero));
  }
  parMap[findCellID][0] += centreWeight;                         // add particleWeights
  parMap[findCellID][1] += pV(radius, scaleVol) * centreWeight;  // add particleVolumes
  parMap[findCellID][2] += pV(radius, scaleVol) * centreWeight;  // add particleVs
}

}  // namespace Foam
