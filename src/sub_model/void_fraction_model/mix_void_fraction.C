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

#include "./divided_void_fraction.h"
#include "./mix_void_fraction.h"

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
}

//! \brief Destructor
mixVoidFraction::~mixVoidFraction() {}

//! \brief 输出空隙率相关信息
void mixVoidFraction::printVoidFractionInfo() const {
  base::MPI_Barrier();
  int nProcs = 0, id = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  for (int i = 0; i < nProcs; ++i) {
    if (id == i) {
      Pout << typeName().c_str() << ":" << endl;
      for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
        Pout << "  findCellIDs[" << index << "]: " << cloud_.findCellIDs()[index] << endl;
        Pout << "  particleOverMeshNumber[" << index << "]: " << cloud_.particleOverMeshNumber()[index] << endl;
        Pout << "  cellIDs[" << index << "]: ";
        for (int j = 0; j < cloud_.particleOverMeshNumber()[index]; ++j) {
          Pout << cloud_.cellIDs()[index][j] << ", ";
        }
        Pout << endl;
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
    base::MPI_Barrier(0.2);
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
  cloud_.pCloud().particleWeights().clear();
  cloud_.pCloud().particleVolumes().clear();
  // set voidFraction
  std::vector<std::unordered_map<int, Foam::vector>> parMaps(cloud_.numberOfParticles());
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkCoarseParticle(index)) {
      // coarse particle
    } else {
      // find and middle particle
      int findCellID = cloud_.findCellIDs()[index];
      if (findCellID >= 0) {
        setVoidFractionForSingleParticle(index, findCellID, parMaps[index]);
      }
    }
  }  // End loop of all particles
  voidFractionNext_.correctBoundaryConditions();
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkCoarseParticle(index)) {
    } else {
      const auto& umap = parMaps[index];
      int meshNumber = umap.size();
      if (meshNumber > 0) {
        // realloc memory according to meshNumber
        cloud_.particleOverMeshNumber()[index] = meshNumber;
        cloud_.pCloud().cellIDs().emplace_back(base::makeShape1(meshNumber), -1);
        cloud_.pCloud().voidFractions().emplace_back(base::makeShape1(meshNumber), 0.0);
        cloud_.pCloud().particleWeights().emplace_back(base::makeShape1(meshNumber), 0.0);
        cloud_.pCloud().particleVolumes().emplace_back(base::makeShape1(meshNumber), 0.0);
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
      } else {
        cloud_.particleOverMeshNumber()[index] = 0;
        cloud_.pCloud().cellIDs().emplace_back();
        cloud_.pCloud().voidFractions().emplace_back();
        cloud_.pCloud().particleWeights().emplace_back();
        cloud_.pCloud().particleVolumes().emplace_back();
      }
    }
  }
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

// //! \brief 设置索引为 index 的单个颗粒的空隙率
// //! \note used for middle particle
// void mixVoidFraction::setVoidFractionForSingleMiddleParticle(const int index, const int findMpiCellID,
//                                                              std::unordered_map<int, Foam::vector>& parMap) {
//   parMap.clear();
//   double radius = cloud_.getRadius(index);                  // 颗粒半径
//   double scaleVol = weight();                               // 颗粒体积因子
//   double scaleRadius = cbrt(porosity());                    // 颗粒半径因子
//   double particleV = pV(radius, scaleVol);                  // 颗粒体积
//   Foam::vector particleCentre = cloud_.getPosition(index);  // 颗粒中心坐标
//   Foam::vector cellPos = Foam::vector::zero;                // 网格坐标
//   radius *= scaleRadius;
//   std::unordered_set<int> set;
//   buildExpandedCellSet(set, findMpiCellID, particleCentre, radius, 1.0);
//   // 获取颗粒覆盖网格的个数
//   int meshNumber = set.size();
//   cloud_.particleOverMeshNumber()[index] = meshNumber;
//   for (int cellID : set) {
//     cellPos = cloud_.mesh().C()[cellID];
//     double core = GaussCore(particleCentre, cellPos, radius);
//     double newAlpha = voidFractionNext_[cellID] - core * particleV;
//     if (newAlpha > alphaMin_) {
//       voidFractionNext_[cellID] = newAlpha;
//     } else {
//       // 如果空隙率低于最低空隙率, 则直接赋值为 alphaMin_
//       voidFractionNext_[cellID] = alphaMin_;
//     }
//     parMap.insert(std::make_pair(cellID, Foam::vector::zero));
//     parMap[cellID][0] += core;             // add particleWeights
//     parMap[cellID][1] = particleV * core;  // add particleVolumes
//     parMap[cellID][2] = particleV * core;  // add particleVs
//   }
// }

}  // namespace Foam
