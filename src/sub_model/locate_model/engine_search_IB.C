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

#include "./engine_search_IB.h"
#include "mathematicalConstants.H"

namespace Foam {

cfdemDefineTypeName(engineSearchIB)

    cfdemCreateNewFunctionAdder(locateModel, engineSearchIB)

    //! \brief Constructor
    engineSearchIB::engineSearchIB(cfdemCloud& cloud, const std::string& derivedTypeName)
    : engineSearch(cloud, derivedTypeName), coef_(2.0) {
  // read properties from dictionary
  verbose_ = subPropsDict_.lookupOrDefault<bool>("verbose", false);
  zSplit_ = subPropsDict_.lookupOrDefault<int>("zSplit", 8);
  xySplit_ = subPropsDict_.lookupOrDefault<int>("xySplit", 16);
  if (zSplit_ < 2 || xySplit_ < 1) {
    FatalError << "Erorr: (zSplit_ < 2 || xySplit_ < 1) with zSplit_ = " << zSplit_ << ", and xySplit_ = " << xySplit_
               << abort(FatalError);
  }
  Info << "xySplit: " << xySplit_ << ", zSplit_: " << zSplit_ << endl;

  // 计算所有 satellite point 个数
  // z 方向: 180 / zSplit_，即 zSplit_ - 1 层
  // xy 方向：360 / xySplit_，即 xySplit_ 层
  // top and bottom: 2
  numberOfSatellitePoints_ = (zSplit_ - 1) * xySplit_ + 2;
  for (int i = 0; i < numberOfSatellitePoints_; ++i) {
    satellitePoints_.push_back(generateSatellitePoint(i));
  }
  base::MPI_Barrier(0.5);

  // set boundBox
  boundBoxPtr_.reset(new boundBox(cloud_.mesh().points(), false));
  if (verbose_) {
    Pout << "Min bounds (x, y, z): " << boundBoxPtr_().min() << endl;
    Pout << "Max bounds (x, y, z): " << boundBoxPtr_().max() << endl;
  }
}

//! \brief Destructor
engineSearchIB::~engineSearchIB() {}

/*!
 * \brief use search engine to get cell id of particle center
 * \param findCellIDs 颗粒覆盖网格的编号
 */
void engineSearchIB::findCell(const base::CITensor1& findCellIDs) const {
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // reset findCellIDs
    findCellIDs[index] = -1;
    // 颗粒中心坐标
    Foam::vector particleCenterPos = cloud_.getPosition(index);
    // 颗粒半径
    double radius = cloud_.getRadius(index);
    Foam::vector periodicPos = particleCenterPos;
    // 判断颗粒中心是否在求解域中
    bool isInside = isInsideRectangularDomain(periodicPos, coef_ * radius);
    if (!isInside && cloud_.checkPeriodicCells()) {
      FatalError << "Error: not support periodic check!" << abort(FatalError);
      // int rangeX = static_cast<int>(cloud_.periodicCheckRange()[0]);
      // int rangeY = static_cast<int>(cloud_.periodicCheckRange()[1]);
      // int rangeZ = static_cast<int>(cloud_.periodicCheckRange()[2]);
      // Foam::vector bMax = globalBox.max();
      // Foam::vector bMin = globalBox.min();
      // for (int xDir = -rangeX; xDir <= rangeX; ++xDir) {
      //   periodicPos[0] = particleCenterPos[0] + xDir * (bMax[0] - bMin[0]);
      //   for (int yDir = -rangeY; yDir <= rangeY; ++yDir) {
      //     periodicPos[1] = particleCenterPos[1] + yDir * (bMax[1] - bMin[1]);
      //     for (int zDir = -rangeZ; zDir <= rangeZ; ++zDir) {
      //       periodicPos[2] = particleCenterPos[2] + zDir * (bMax[2] - bMin[2]);
      //       isInside = isInsideRectangularDomain(periodicPos, coef_ * radius);
      //       if (isInside) { break; }
      //     }
      //     if (isInside) { break; }
      //   }
      //   if (isInside) { break; }
      // } // end of check periodic
    }
    if (isInside) {
      // findSingleCell is defined in engineSearch and return the found cell id.
      findCellIDs[index] = findSingleCell(periodicPos, -1);
      // not found
      if (findCellIDs[index] < 0) {
        label altStartID = -1;
        // 遍历当前颗粒的所有 satellitePoint
        for (size_t i = 0; i < satellitePoints_.size(); ++i) {
          Foam::vector pos = getSatellitePointPos(index, static_cast<int>(i));
          isInside = isInsideRectangularDomain(pos, Foam::SMALL);
          // 如果当前 satellite point is in rectangular domain
          if (isInside) {
            altStartID = findSingleCell(pos, -1);
          }
          if (cloud_.checkPeriodicCells()) {
            FatalError << "Error: not support periodic check!" << abort(FatalError);
          }
          if (altStartID >= 0) {
            // found cell id
            findCellIDs[index] = altStartID;
            break;
          }
        }
      }
    } else {
      findCellIDs[index] = -1;
      Pout << "Found particle " << index
           << " not in the CFD domian, this could means that the particle and fluid has no coupling." << endl;
    }
  }  // end loop of particle
}

/*!
 * \brief 判断 pos 是否位于求解域中
 * \param pos 颗粒中心位置
 * \param offsetValue 位置偏移量
 */
bool engineSearchIB::isInsideRectangularDomain(const Foam::vector& pos, double offsetValue) const {
  Foam::vector offset(offsetValue, offsetValue, offsetValue);
  boundBox box(boundBoxPtr_().min() - offset, boundBoxPtr_().max() + offset);
  return box.contains(pos);
}

/*!
 * \brief generate satellite point according to index
 * \param index satellite point index
 */
Foam::vector engineSearchIB::generateSatellitePoint(int index) const {
  double degToRad = M_PI / 180.0;
  double thetaSize = 180.0 / zSplit_;
  double phiSize = 360.0 / xySplit_;
  Foam::vector res(0.0, 0.0, 0.0);
  // one satellite point at bottom, one satellite point at top
  if (index == 0) {
    res[2] = 1.0;
  } else if (index == 1) {
    res[2] = -1.0;
  } else {
    int thetaIndex = (index - 2) / xySplit_;
    int phiIndex = (index - 2) % xySplit_;
    double theta = degToRad * thetaSize * (thetaIndex + 1);
    double phi = degToRad * phiSize * phiIndex;
    res[0] = sin(theta) * cos(phi);
    res[1] = sin(theta) * sin(phi);
    res[2] = cos(theta);
  }
  return res;
}

/*!
 * \brief get satellite point position
 * \param index particle index
 * \param satellitePointIndex satellite point index
 */
inline Foam::vector engineSearchIB::getSatellitePointPos(int index, int satellitePointIndex) const {
  double radius = cloud_.getRadius(index);
  Foam::vector particleCenterPos = cloud_.getPosition(index);
  return particleCenterPos + radius * satellitePoints_[satellitePointIndex];
}

}  // namespace Foam
