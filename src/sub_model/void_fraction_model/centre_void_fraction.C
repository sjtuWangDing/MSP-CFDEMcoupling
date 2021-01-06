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

#include "./centre_void_fraction.h"

namespace Foam {

cfdemDefineTypeName(centreVoidFraction);

cfdemCreateNewFunctionAdder(voidFractionModel, centreVoidFraction);

//! \brief Constructor
centreVoidFraction::centreVoidFraction(cfdemCloud& cloud) : voidFractionModel(cloud), alphaMin_(0.0) {
  weight_ = subPropsDict_.lookupOrDefault<double>("weight", 1.0);
  alphaMin_ = subPropsDict_.lookupOrDefault<double>("alphaMin", 0.0);
}

//! \brief Destructor
centreVoidFraction::~centreVoidFraction() {}

void centreVoidFraction::setVoidFraction() {
  // reset field
  resetVoidFraction();
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
    setVoidFractionForSingleParticle(index, parMaps[index]);
  }  // End loop of all particles
  voidFractionNext_.correctBoundaryConditions();
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    const auto& umap = parMaps[index];
    int meshNumber = umap.size();
    CHECK(meshNumber == 0 || meshNumber == 1) << ": Wrong setting void fraction for particle with index: " << index;
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

//! \brief 设置索引为 index 的单个颗粒的空隙率
void centreVoidFraction::setVoidFractionForSingleParticle(const int index,
                                                          std::unordered_map<int, Foam::vector>& parMap) {
  // 获取颗粒半径
  double radius = cloud_.getRadius(index);
  // 获取到在当前 processor 上颗粒覆盖的某一个网格编号
  int findCellID = cloud_.findCellIDs()[index];
  // 使用 scaleVol 作为颗粒体积因子
  double scaleVol = weight();
  // particle centre is in domain
  if (findCellID >= 0) {
    double findCellVol = cloud_.mesh().V()[findCellID];
    double particleVol = pV(radius, scaleVol);
    // 更新空隙率场
    double newAlpha = voidFractionNext_[findCellID] - particleVol / findCellVol;
    voidFractionNext_[findCellID] = newAlpha > alphaMin_ ? newAlpha : alphaMin_;
    parMap.insert(std::make_pair(findCellID, Foam::vector(1.0, particleVol, particleVol)));
  }
}

}  // namespace Foam
