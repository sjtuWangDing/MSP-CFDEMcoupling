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
  Foam::ParticleCloud
\*---------------------------------------------------------------------------*/

#ifndef __PARTICLE_CLOUD_H__
#define __PARTICLE_CLOUD_H__

#include "base/tensor/tensor.h"
#include "cloud/cfdem_base.h"

namespace Foam {

class cfdemCloud;

//! \brief particle dimension type enum
enum ParType { kFine = 0, kMiddle, kCoarse };

//! \brief The basic class for particles
class ParticleCloud {
  friend class cfdemCloud;

 public:
  ParticleCloud(int numberOfParticles = -1)
      : numberOfParticles_(numberOfParticles),
        numberOfParticlesChanged_(false),
        radiiPtr_(nullptr),
        cdsPtr_(nullptr),
        positionsPtr_(nullptr),
        velocitiesPtr_(nullptr),
        initVelocitiesPtr_(nullptr),
        angularVelocitiesPtr_(nullptr),
        DEMForcesPtr_(nullptr),
        DEMTorquesPtr_(nullptr),
        fluidVelPtr_(nullptr) {}

  inline int& numberOfParticles() { return numberOfParticles_; }

  inline bool& numberOfParticlesChanged() { return numberOfParticlesChanged_; }

  inline base::CITensor1& particleOverMeshNumber() { return particleOverMeshNumber_; }

  inline base::CITensor1& findCellIDs() { return findCellIDs_; }

  inline base::CDTensor1& dimensionRatios() { return dimensionRatios_; }

  inline base::CDTensor2& fAcc() { return fAcc_; }

  inline base::CDTensor2& impForces() { return impForces_; }

  inline base::CDTensor2& expForces() { return expForces_; }

  inline std::vector<base::CITensor1>& cellIDs() { return cellIDs_; }

  inline std::vector<base::CDTensor1>& voidFractions() { return voidFractions_; }

  inline std::vector<base::CDTensor1>& volumeFractions() { return volumeFractions_; }

  inline std::vector<base::CDTensor1>& particleWeights() { return particleWeights_; }

  inline std::vector<base::CDTensor1>& particleVolumes() { return particleVolumes_; }

  inline double**& radiiPtr() { return radiiPtr_; }

  inline double**& cdsPtr() { return cdsPtr_; }

  inline double**& positionsPtr() { return positionsPtr_; }

  inline double**& velocitiesPtr() { return velocitiesPtr_; }

  inline double**& initVelocitiesPtr() { return initVelocitiesPtr_; }

  inline double**& angularVelocitiesPtr() { return angularVelocitiesPtr_; }

  inline double**& DEMForcesPtr() { return DEMForcesPtr_; }

  inline double**& DEMTorquesPtr() { return DEMTorquesPtr_; }

  inline double**& fluidVelPtr() { return fluidVelPtr_; }

  inline base::CDExTensor1& cds() { return cds_; }

  inline base::CDExTensor1& radii() { return radii_; }

  inline base::CDExTensor2& positions() { return positions_; }

  inline base::CDExTensor2& velocities() { return velocities_; }

  inline base::CDExTensor2& initVelocities() { return initVelocities_; }

  inline base::CDExTensor2& angularVelocities() { return angularVelocities_; }

  inline base::CDExTensor2& DEMForces() { return DEMForces_; }

  inline base::CDExTensor2& DEMTorques() { return DEMTorques_; }

  inline base::CDExTensor2& fluidVel() { return fluidVel_; }

  inline double getRadius(int index) const {
    CHECK_EQ(numberOfParticles_, static_cast<int>(radii_.size(0))) << __func__ << ": Number of particle is not match";
    CHECK_GT(radii_[index], Foam::SMALL) << __func__ << ": Radius of particle " << index << " is < "
                                         << "Foam::SMALL";
    return radii_[index];
  }

  inline double getDimensionRatio(int index) const {
    CHECK_EQ(numberOfParticles_, static_cast<int>(dimensionRatios_.size(0)))
        << __func__ << ": Number of particle is not match";
    return dimensionRatios_[index];
  }

  inline Foam::vector getPosition(int index) const {
    CHECK_EQ(numberOfParticles_, static_cast<int>(positions_.size(0)))
        << "getPosition: Number of particle is not match";
    Foam::vector pos = Foam::vector::zero;
#pragma unroll
    for (int i = 0; i < 3; ++i) {
      pos[i] = positions_[index][i];
    }
    return pos;
  }

  inline Foam::vector getVelocity(int index) const {
    CHECK_EQ(numberOfParticles_, static_cast<int>(velocities_.size(0)))
        << "getVelocity: Number of particle is not match";
    Foam::vector vel = Foam::vector::zero;
#pragma unroll
    for (int i = 0; i < 3; ++i) {
      vel[i] = velocities_[index][i];
    }
    return vel;
  }

  inline Foam::vector getInitVelocity(int index) const {
    CHECK_EQ(numberOfParticles_, static_cast<int>(initVelocities_.size(0)))
        << "getInitVelocity: Number of particle is not match";
    Foam::vector vel = Foam::vector::zero;
#pragma unroll
    for (int i = 0; i < 3; ++i) {
      vel[i] = initVelocities_[index][i];
    }
    return vel;
  }

  inline Foam::vector getAngularVelocity(int index) const {
    CHECK_EQ(numberOfParticles_, static_cast<int>(angularVelocities_.size(0)))
        << "getAngularVelocity: Number of particle is not match";
    Foam::vector vel = Foam::vector::zero;
#pragma unroll
    for (int i = 0; i < 3; ++i) {
      vel[i] = angularVelocities_[index][i];
    }
    return vel;
  }

  inline void setNumberOfParticles(int number) {
    CHECK_GE(number, 0) << "setNumberOfParticles: number of particle is <= 0";
    numberOfParticlesChanged_ = numberOfParticles_ != number ? true : false;
    numberOfParticles_ = number;
  }

 private:
  //! \brief 颗粒总数
  int numberOfParticles_;
  //! \brief 颗粒总数是否发生变化
  bool numberOfParticlesChanged_;
#if CFDEM_USE_TENSOR
  base::CITensor1 particleOverMeshNumber_;
  base::CITensor1 findCellIDs_;
  base::CDTensor1 dimensionRatios_;
  base::CDTensor2 fAcc_;
  base::CDTensor2 impForces_;
  base::CDTensor2 expForces_;
  std::vector<base::CITensor1> cellIDs_;
  std::vector<base::CDTensor1> voidFractions_;
  std::vector<base::CDTensor1> volumeFractions_;
  std::vector<base::CDTensor1> particleWeights_;
  std::vector<base::CDTensor1> particleVolumes_;
  double** radiiPtr_;
  double** cdsPtr_;
  double** positionsPtr_;
  double** velocitiesPtr_;
  double** initVelocitiesPtr_;
  double** angularVelocitiesPtr_;
  double** DEMForcesPtr_;
  double** DEMTorquesPtr_;
  double** fluidVelPtr_;
  base::CDExTensor1 radii_;
  base::CDExTensor1 cds_;
  base::CDExTensor2 positions_;
  base::CDExTensor2 velocities_;
  base::CDExTensor2 initVelocities_;
  base::CDExTensor2 angularVelocities_;
  base::CDExTensor2 DEMForces_;
  base::CDExTensor2 DEMTorques_;
  base::CDExTensor2 fluidVel_;
#else
  //! \brief 颗粒半径
  //! \note radii_[index][0]
  double** radii_;

  //! \brief 颗粒体积
  //! \note volumes_[index][0]
  double** volumes_;

  //! \brief 颗粒阻力系数
  //         当 forceSubModel implForceDEM 为 true 时候，该 cds_ 和 fluidVel_ 传递给 DEM 求解器，
  //         DEM 求解器使用颗粒中心处的流体速度与阻力系数一起计算颗粒受到的阻力
  //! \note cds_[index][0]
  double** cds_;

  //! \brief 颗粒位置矢量
  //! \note positions_[index][0 ~ 2]
  double** positions_;

  //! \brief 颗粒速度
  //! \note velocities_[index][0 ~ 2]
  double** velocities_;

  //! \brief 小颗粒的中心处的流体速度
  //         当 forceSubModel implForceDEM 为 true 时候，该 cds_ 和 fluidVel_ 传递给 DEM 求解器，
  //         DEM 求解器使用颗粒中心处的流体速度与阻力系数一起计算颗粒受到的阻力
  //! \note fluidVel_[index][0 ~ 2]
  double** fluidVel_;

  //! \brief 在一个耦合时间间隔中，小颗粒所受的阻力的累计，将在下个耦合步中，传递给流场
  //! \note fAcc_[index][0 ~ 2] ??????????
  double** fAcc_;

  //! \brief 颗粒对流体的隐式作用力
  //! \note impForces_[index][0 ~ 2]
  double** impForces_;

  //! \brief 颗粒对流体的显式作用力
  //! \note expForces_[index][0 ~ 2]
  double** expForces_;

  //! \brief 流体对颗粒的总作用力，如果不设置 cds_，则
  //! \note DEMForces_[index][0 ~ 2]
  double** DEMForces_;

  //! \brief 空隙率
  //! \note voidFractions_[index][subcell]
  double** voidFractions_;

  //! \brief 颗粒覆盖的所有网格的编号
  //! \note cellIDs_[index][subcell]
  double** cellIDs_;

  //! \brief 颗粒对所覆盖网格的影响系数
  //         如果使用 divided 空隙率模型，则对颗粒覆盖的某个网格 subcell，
  //         只要有一个颗粒标志点在网格中，particleWeights_[index][subcell] += 1.0 / 29.0
  //! \note particleWeights_[index][subcell]
  double** particleWeights_;

  //! \brief 颗粒对所覆盖网格的覆盖体积
  //         如果使用 divided 空隙率模型，则对颗粒覆盖的某个网格 subcell，
  //         只要有一个颗粒标志点在网格中，particleVolumes_[index][subcell] += (1.0 / 29.0) * 颗粒体积
  //! \note particleVolumes_[index][subcell]
  double** particleVolumes_;
#endif  // CFDEM_CLOUD_USE_TENSOR
};

}  // namespace Foam

#endif  // __PARTICLE_CLOUD_H__
