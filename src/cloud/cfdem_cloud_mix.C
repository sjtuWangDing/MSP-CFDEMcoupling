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

#include "cloud/cfdem_cloud_mix.h"
#include "sub_model/averaging_model/averaging_model.h"
#include "sub_model/force_model/force_model.h"
#include "sub_model/locate_model/locate_model.h"
#include "sub_model/mom_couple_model/implicit_couple.h"
#include "sub_model/mom_couple_model/mom_couple_model.h"
#include "sub_model/void_fraction_model/void_fraction_model.h"

namespace Foam {

cfdemDefineTypeName(cfdemCloudMix);

//! \brief Constructed from mesh
cfdemCloudMix::cfdemCloudMix(const fvMesh& mesh) : cfdemCloud(mesh) {}

//! \brief Destructor
cfdemCloudMix::~cfdemCloudMix() {}

//! \brief 重新分配内存
void cfdemCloudMix::reallocate() {
  int number = numberOfParticles();
  // allocate memory of data exchanged with liggghts
  dataExchangeM().realloc(parCloud_.radii(), base::makeShape1(number), parCloud_.radiiPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.cds(), base::makeShape1(number), parCloud_.cdsPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.positions(), base::makeShape2(number, 3), parCloud_.positionsPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.velocities(), base::makeShape2(number, 3), parCloud_.velocitiesPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.DEMForces(), base::makeShape2(number, 3), parCloud_.DEMForcesPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.fluidVel(), base::makeShape2(number, 3), parCloud_.fluidVelPtr(), 0.0);
  // allocate memory of data not exchanged with liggghts
  parCloud_.particleOverMeshNumber() = std::move(base::CITensor1(base::makeShape1(number), 0));
  parCloud_.dimensionRatios() = std::move(base::CDTensor1(base::makeShape1(number), -1.0));
  parCloud_.impForces() = std::move(base::CDTensor2(base::makeShape2(number, 3), 0.0));
  parCloud_.expForces() = std::move(base::CDTensor2(base::makeShape2(number, 3), 0.0));
  parCloud_.particleRootProcIDs() = std::move(base::CITensor1(base::makeShape1(number), -1));
  parCloud_.findCellIDs() = std::move(base::CITensor1(base::makeShape1(number), -1));
  parCloud_.findMpiCellIDs() = std::move(base::CITensor1(base::makeShape1(number), -1));
  parCloud_.findExpandedCellIDs() = std::move(base::CITensor1(base::makeShape1(number), -1));
}

void cfdemCloudMix::printParticleInfo() const {
  int nProcs = 0, id = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  base::MPI_Barrier();
  if (0 == id) {
    for (int index = 0; index < numberOfParticles(); ++index) {
      Pout << "  position[" << index << "]: " << positions()[index][0] << ", " << positions()[index][1] << ", "
           << positions()[index][2] << endl;
    }
    for (int index = 0; index < numberOfParticles(); ++index) {
      Pout << "  velocity[" << index << "]: " << velocities()[index][0] << ", " << velocities()[index][1] << ", "
           << velocities()[index][2] << endl;
    }
  }
  base::MPI_Barrier();
  for (int index = 0; index < numberOfParticles(); ++index) {
    Pout << "  dimensionRatio[" << index << "]: " << dimensionRatios()[index] << endl;
  }
  base::MPI_Barrier();
  for (int index = 0; index < numberOfParticles(); ++index) {
    Pout << "  DEMForce[" << index << "]: " << DEMForces()[index][0] << ", " << DEMForces()[index][1] << ", "
         << DEMForces()[index][2] << endl;
  }
  base::MPI_Barrier();
  for (int index = 0; index < numberOfParticles(); ++index) {
    Pout << "  impForce[" << index << "]: " << impForces()[index][0] << ", " << impForces()[index][1] << ", "
         << impForces()[index][2] << endl;
  }
  base::MPI_Barrier();
  voidFractionM().printVoidFractionInfo();
}

/*!
 * \brief 更新函数
 * \note used for cfdemSolverPiso
 * \param U      <[in] 流体速度场
 * \param voidF  <[in, out] 小颗粒空隙率场
 * \param Us     <[in, out] 局部平均小颗粒速度场
 * \param Ksl    <[in, out] 动量交换场
 */
void cfdemCloudMix::evolve(volVectorField& U, volScalarField& voidF, volVectorField& Us, volScalarField& Ksl) {
  Info << "/ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /" << endl;
  // 检查当前流体时间步是否同时也是耦合时间步
  validCouplingStep_ = dataExchangeM().checkValidCouplingStep();
  if (validCouplingStep_) {
    // 创建用于记录 coupling time step counter
    auto pCounter = std::make_shared<dataExchangeModel::CouplingStepCounter>(dataExchangeM());
    // couple(): run liggghts command and get number of particle
    setNumberOfParticles(dataExchangeM().couple());
    // reset field
    resetField();
    // realloc memory
    reallocate();
    // 获取 DEM data
    getDEMData();
    // 获取颗粒中心所在的网格索引
    locateM().findCell(parCloud_.findCellIDs());
    // 获取颗粒覆盖当前处理器上的某一个网格索引，如果索引为 -1，则表示颗粒不覆盖当前处理器
    locateM().findMpiCell(parCloud_.findMpiCellIDs());
    // 计算颗粒尺度
    voidFractionM().getDimensionRatiosForMix(parCloud_.dimensionRatios());
    // 计算扩展颗粒覆盖当前处理器上的某一个网格索引
    // 必须位于计算颗粒尺度之后，因为需要判断颗粒是否为 middle
    locateM().findExpandedCell(parCloud_.findExpandedCellIDs(), expandedCellScale());
    // 计算颗粒空隙率
    voidFractionM().setVoidFraction();
    voidF = voidFractionM().voidFractionInterp();
    // 计算局部平局颗粒速度场
    averagingM().setVectorFieldAverage(averagingM().UsNext(), averagingM().UsWeightField(), velocities(),
                                       particleWeights());
    Us = averagingM().UsInterp();
    // global force init
    globalF().initBeforeSetForce();
    // 计算流体对颗粒的作用力
    for (const auto& ptr : forceModels_) {
      ptr->setForce();
    }
    // global force end
    globalF().endAfterSetForce();
    // 计算局部累加的流体作用力场
    averagingM().setVectorFieldSum(globalF().impParticleForce(), impForces(), particleWeights());
    // write DEM data
    giveDEMData();
    // get shared ptr of implicitCouple model and update Ksl field
    std::shared_ptr<momCoupleModel> sPtr = momCoupleModels_[implicitCouple::cTypeName()];
    Ksl = sPtr->impMomSource();
    Ksl.correctBoundaryConditions();
    printParticleInfo();
  }
  Info << __func__ << " - done\n" << endl;
}

}  // namespace Foam
