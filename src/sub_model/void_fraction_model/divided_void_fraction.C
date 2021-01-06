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

#include <unordered_map>
#include <vector>
#include "./divided_void_fraction.h"
#include "mpi.h"

namespace Foam {

cfdemDefineTypeName(dividedVoidFraction);

cfdemCreateNewFunctionAdder(voidFractionModel, dividedVoidFraction);

//! \brief Constructor
dividedVoidFraction::dividedVoidFraction(cfdemCloud& cloud)
    : voidFractionModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      verbose_(false),
      alphaMin_(0.0),
      tooMuch_(0.0) {
  weight_ = subPropsDict_.lookupOrDefault<double>("weight", 1.0);
  porosity_ = subPropsDict_.lookupOrDefault<double>("porosity", 1.0);
  verbose_ = subPropsDict_.lookupOrDefault<bool>("verbose", false);
  alphaMin_ = subPropsDict_.lookupOrDefault<double>("alphaMin", 0.0);
  // 计算标志点相对颗粒中心的坐标
  [this](void) -> void {
    int idx = 0;
    // 设置中心标志点偏移
    for (int i = 0; i < 3; ++i) {
      offsets_[idx][i] = 0.0;
    }
    idx += 1;
    // 计算两个半径, 这两个半径构成的球面将颗粒分为 1 : 14 : 14 这三个部分, 这三个部分分别是: 中心球体a, 球面层 b,
    // 球面层 c
    // 第一个部分是一个中心与粒子中心重合的球体。该子颗粒的半径 r1 计算如下: V(r1) / V(R) = r1 ^ 3 / R ^ 3 = 1 / 29 ==>
    // r1 = R * (1/29)^(1/3)
    // 其余体积是一个球面层, 必须分为两个等体积球面层。这两个球面层之间的径向边界位置很容易获得: V(r1) / V(R) = r1 ^ 3 /
    // R ^ 3 = 15 / 29 ==> r1 = R * (15/29)^(1/3)
    double r1 = cbrt(1.0 / double(numberOfMarkerPoints_));
    double r2 = cbrt(15.0 / double(numberOfMarkerPoints_));

    // 将球面层 b 和 c 分别等体积分成 14 个部分
    // 第一球面层中各体积径向质心点的位置是: r(s,1) = 0.62761 * R
    // 第二球面层中各体积径向质心点的位置是: r(s,2) = 0.90853 * R
    scalar r[2] = {0.75 * (r2 * r2 * r2 * r2 - r1 * r1 * r1 * r1) / (r2 * r2 * r2 - r1 * r1 * r1),
                   0.75 * (1.0 - r2 * r2 * r2 * r2) / (1.0 - r2 * r2 * r2)};

    for (int ir = 0; ir < 2; ++ir) {
      // 通过球坐标系找到 8 个标志点
      for (scalar zeta = 0.25 * M_PI; zeta < 2.0 * M_PI; zeta += 0.5 * M_PI) {
        for (scalar theta = 0.25 * M_PI; theta < M_PI; theta += 0.5 * M_PI) {
          offsets_[idx][0] = r[ir] * Foam::sin(theta) * Foam::cos(zeta);
          offsets_[idx][1] = r[ir] * Foam::sin(theta) * Foam::sin(zeta);
          offsets_[idx][2] = r[ir] * Foam::cos(theta);
          idx += 1;
        }
      }
      // 通过笛卡尔坐标系找到 6 个标志点(x y z 方向各有两个)
      for (int j = -1; j <= 1; j += 2) {
        offsets_[idx][0] = r[ir] * static_cast<double>(j);
        offsets_[idx][1] = 0.0;
        offsets_[idx][2] = 0.0;
        idx += 1;

        offsets_[idx][0] = 0.0;
        offsets_[idx][1] = r[ir] * static_cast<double>(j);
        offsets_[idx][2] = 0.0;
        idx += 1;

        offsets_[idx][0] = 0.0;
        offsets_[idx][1] = 0.0;
        offsets_[idx][2] = r[ir] * static_cast<double>(j);
        idx += 1;
      }
    }
  }();
}

//! \brief Destructor
dividedVoidFraction::~dividedVoidFraction() {}

//! \brief 输出空隙率相关信息
void dividedVoidFraction::printVoidFractionInfo() const {
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

void dividedVoidFraction::setVoidFraction() {
  // reset field
  resetVoidFraction();
  // reset particleOverMeshNumber
  base::fillTensor(cloud_.particleOverMeshNumber(), 0);
  // clear before set
  cloud_.pCloud().cellIDs().clear();
  cloud_.pCloud().voidFractions().clear();
  cloud_.pCloud().particleWeights().clear();
  cloud_.pCloud().particleVolumes().clear();
  // reset tooMuch_
  tooMuch_ = 0.0;
  // set voidFraction
  std::vector<std::unordered_map<int, Foam::vector>> parMaps(cloud_.numberOfParticles());
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    setVoidFractionForSingleParticle(index, parMaps[index]);
  }  // End loop of all particles
  voidFractionNext_.correctBoundaryConditions();
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    const auto& umap = parMaps[index];
    int meshNumber = umap.size();
    CHECK(meshNumber >= 0 && meshNumber <= numberOfMarkerPoints_)
        << ": Wrong setting void fraction for particle with index: " << index;
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

void dividedVoidFraction::setVoidFractionForSingleParticle(const int index,
                                                           std::unordered_map<int, Foam::vector>& parMap) {
  parMap.clear();
  // 获取颗粒半径
  double radius = cloud_.getRadius(index);
  // 获取颗粒中心坐标
  Foam::vector particleCentre = cloud_.getPosition(index);
  // 获取到在当前 processor 上颗粒覆盖的某一个网格编号
  int findCellID = cloud_.findCellIDs()[index];
  // 使用 scaleVol 作为颗粒体积因子
  double scaleVol = weight();
  double scaleRadius = cbrt(porosity());
  // 获取颗粒体积的 1 / 29
  double particleVolume = pV(radius, scaleVol) / static_cast<double>(numberOfMarkerPoints_);
  // 先计算体积, 然后再用 scaleRadius 乘以半径, 以保证在计算空隙率时候颗粒尺寸不变
  radius *= scaleRadius;
  // particle centre is in domain
  if (findCellID >= 0) {
    // 统计有效标记点的个数
    int findSubCellNum = 0;
    // 遍历所有标志点，除了颗粒中心处的标志点
    for (int i = 1; i < numberOfMarkerPoints_; ++i) {
      // 获取标志点的绝对位置
      Foam::vector subPosition = particleCentre + radius * offsets_[i];
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
        double newAlpha = voidFractionNext_[subCellId] - particleVolume / subCellVol;
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
        parMap[subCellId][0] += 1.0 / static_cast<double>(numberOfMarkerPoints_);  // add particleWeights
        parMap[subCellId][1] += particleVolume;                                    // add particleVolumes
        parMap[subCellId][2] += particleVolume;                                    // add particleVs
      }
    }  // End of marker point loop
    // 处理颗粒中心处的标志点, 将丢失的标志点对空隙率的影响作用在颗粒中心处的标志点上
    // set source for particle center
    double findCellVol = cloud_.mesh().V()[findCellID];
    double centreWeight = (numberOfMarkerPoints_ - findSubCellNum) * (1.0 / numberOfMarkerPoints_);
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
}

}  // namespace Foam
