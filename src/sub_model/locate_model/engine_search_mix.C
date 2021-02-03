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

#include "./engine_search_mix.h"
#include "mpi.h"

namespace Foam {

cfdemDefineTypeName(engineSearchMix);

cfdemCreateNewFunctionAdder(locateModel, engineSearchMix);

//! \brief Constructor
engineSearchMix::engineSearchMix(cfdemCloud& cloud, const std::string& derivedTypeName)
    : engineSearchIB(cloud, derivedTypeName) {}

//! \brief Destructor
engineSearchMix::~engineSearchMix() {}

/*!
 * \brief use search engine to get cell id of particle center
 * \param findCellIDs 颗粒覆盖网格的编号
 */
void engineSearchMix::findCell(const base::CITensor1& findCellIDs) const {
  // 如果在当前耦合时间步中网格更新或者模型允许耦合时间步长 > 流体时间步长，那么就可能在非耦合的流体时间步中更新网格
  // 所以如果 allowUseSubCFDTimeStep 为 true，那么保险起见，在每个耦合时间步中都修正 searchEngine
  // 只需要在三个搜索函数中调用一次 searchEngine_.correct() 即可
  if (cloud_.meshHasUpdated() || cloud_.allowUseSubCFDTimeStep()) {
    dynamic_cast<engineSearchIB*>(const_cast<engineSearchMix*>(this))->searchEngine().correct();
  }
  engineSearch::findCell(findCellIDs);
  base::MPI_Info("engineSearchMix: findCell - done", verbose_);
}

/*!
 * \brief use search engine to get id of cell which covered by processor
 * \param findMpiCellIDs 颗粒覆盖当前求解器上任意一个网格的编号
 * \param scale 颗粒半径尺度系数
 */
void engineSearchMix::findMpiCell(const base::CITensor1& findMpiCellIDs) const {
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    findMpiCellKernel(index, findMpiCellIDs);
    if (-1 == findMpiCellIDs[index] && 0 == index && verbose_) {
      Pout << __func__ << ": Find no mesh covered by particle " << index
           << ", this could means the particle is not in the CFD domian." << endl;
    }
  }
  base::MPI_Info("engineSearchMix: findMpiCell - done", verbose_);
}

/*!
 * \brief use search engine to get id of cell which covered by processor
 * \param findExpandedCellIDs 被扩展的颗粒覆盖当前求解器上任意一个网格的编号
 */
void engineSearchMix::findExpandedCell(const base::CITensor1& findExpandedCellIDs, const double scale) const {
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkMiddleParticle(index)) {
      findMpiCellKernel(index, findExpandedCellIDs, scale);
    }
  }
  base::MPI_Info("engineSearchMix: findExpandedCell - done", verbose_);
}

void engineSearchMix::findMpiCellKernel(const int index, const base::CITensor1& findMpiCellIDs,
                                        const double scale /* = 1.0 */) const {
  double radius = cloud_.getRadius(index);                     // 颗粒半径
  Foam::vector particleCenterPos = cloud_.getPosition(index);  // 颗粒中心坐标
  Foam::vector spPos = Foam::vector::zero;                     // 标志点坐标
  findMpiCellIDs[index] = -1;
  // 判断颗粒中心是否在求解域中
  bool isInside = isInsideRectangularDomain(particleCenterPos, coef_ * radius);
  if (!isInside && cloud_.checkPeriodicCells()) {
    FatalError << "Error: not support periodic check!" << abort(FatalError);
  }
  if (isInside) {
    findMpiCellIDs[index] = findSingleCell(particleCenterPos, -1);
    if (findMpiCellIDs[index] < 0) {  // not found
      int altStartID = -1;
      // 遍历当前颗粒的所有 satellitePoint
      for (size_t i = 0; i < satellitePoints_.size(); ++i) {
        spPos = getSatellitePointPos(index, static_cast<int>(i), scale);
        isInside = isInsideRectangularDomain(spPos, Foam::SMALL);
        // 如果当前 satellite point is in rectangular domain
        if (isInside) {
          altStartID = findSingleCell(spPos, -1);
        }
        if (cloud_.checkPeriodicCells()) {
          FatalError << "Error: not support periodic check!" << abort(FatalError);
        }
        if (altStartID >= 0) {
          // found cell id
          findMpiCellIDs[index] = altStartID;
          break;
        }
      }
    }
  }
}

}  // namespace Foam