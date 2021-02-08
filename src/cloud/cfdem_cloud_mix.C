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
#include "dynamicFvMesh.H"
#include "sub_model/averaging_model/averaging_model.h"
#include "sub_model/force_model/force_model.h"
#include "sub_model/locate_model/locate_model.h"
#include "sub_model/mom_couple_model/implicit_couple.h"
#include "sub_model/mom_couple_model/mom_couple_model.h"
#include "sub_model/void_fraction_model/void_fraction_model.h"

namespace Foam {

cfdemDefineTypeName(cfdemCloudMix);

//! \brief Constructed from mesh
cfdemCloudMix::cfdemCloudMix(const fvMesh& mesh)
    : cfdemCloud(mesh),
      pRefCell_(readLabel(mesh.solutionDict().subDict("PISO").lookup("pRefCell"))),
      pRefValue_(readScalar(mesh.solutionDict().subDict("PISO").lookup("pRefValue"))) {}

//! \brief Destructor
cfdemCloudMix::~cfdemCloudMix() {}

//! \brief 重新分配内存
void cfdemCloudMix::reallocate() {
  int number = numberOfParticles();
  // allocate memory of data exchanged with liggghts
  dataExchangeM().realloc(parCloud_.radii(), base::makeShape1(number), parCloud_.radiiPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.positions(), base::makeShape2(number, 3), parCloud_.positionsPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.velocities(), base::makeShape2(number, 3), parCloud_.velocitiesPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.DEMForces(), base::makeShape2(number, 3), parCloud_.DEMForcesPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.cds(), base::makeShape1(number), parCloud_.cdsPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.fluidVel(), base::makeShape2(number, 3), parCloud_.fluidVelPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.angularVelocities(), base::makeShape2(number, 3), parCloud_.angularVelocitiesPtr(),
                          0.0);
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

void cfdemCloudMix::getDEMData() {
  cfdemCloud::getDEMData();
  dataExchangeM().getData("omega", "vector-atom", parCloud_.angularVelocitiesPtr());
}

void cfdemCloudMix::printParticleInfo() const {
  base::MPI_Barrier();
  auto info = [this](const int index) {
    int rootId = particleRootProcIDs()[index];
    Foam::vector position(positions()[index][0], positions()[index][1], positions()[index][2]);
    Foam::vector velocity(velocities()[index][0], velocities()[index][1], velocities()[index][2]);
    Foam::vector DEMForce(DEMForces()[index][0], DEMForces()[index][1], DEMForces()[index][2]);
    Foam::vector impForce(impForces()[index][0], impForces()[index][1], impForces()[index][2]);
    if (rootId == base::procId()) {
      Pout << "  position[" << index << "]: " << position << endl;
      Pout << "  velocity[" << index << "]: " << velocity << endl;
      Pout << "  dimensionRatio[" << index << "]: " << dimensionRatios()[index] << endl;
      if (mag(DEMForce) > Foam::SMALL) {
        Pout << "  DEMForce[" << index << "]: " << DEMForce << endl;
      }
      if (mag(impForce) > Foam::SMALL) {
        Pout << "  impForce[" << index << "]: " << impForce << endl;
      }
    }
  };
  std::once_flag onceC, onceM, onceF;
  for (int index = 0; index < numberOfParticles(); ++index) {
    if (checkCoarseParticle(index)) {
      std::call_once(onceC, info, index);
    } else if (checkMiddleParticle(index)) {
      std::call_once(onceM, info, index);
    } else {
      std::call_once(onceF, info, index);
    }
  }
  base::MPI_Barrier();
  // voidFractionM().printVoidFractionInfo();
}

/*!
 * \brief 更新网格，如果 mesh 是 Foam::dynamicRefineFvMesh 类型，则更新网格，如果是 Foam::staticFvMesh
 * 或者其他类型，则不更新
 */
void cfdemCloudMix::updateMesh(volScalarField& interface) {
  base::MPI_Info("update mesh...", true);
  fvMesh& mesh = const_cast<fvMesh&>(mesh_);
  try {
    // 采用动态网格
    dynamicFvMesh& dyMesh = dynamic_cast<dynamicFvMesh&>(mesh);
    // 设置 interface，范围是 coarse particle 的半径的 refineMeshSkin 倍
    setInterface(interface, refineMeshSkin());
    interface.correctBoundaryConditions();
    // 如果 dyMesh.update() 返回 true，则表明 mesh 被更新了
    // 在findCell函数中，searchEngine会调用correct()函数修正
    setMeshHasUpdated(dyMesh.update());
    base::MPI_Barrier();
    // 这里必须调用 mesh.C(), 否则在调用 mesh.C()[cellID] 获取网格中心坐标的时候报错：An error occurred in MPI_Waitall
    Pout << "Mesh number in current Proc: " << dyMesh.C().size() << endl;
  } catch (const std::bad_cast& ex) {
    Info << "Not use dynamicRefineFvMesh, no need to update mesh" << endl;
  }
  base::MPI_Info("update mesh - done", true);
  return;
}

//! \brief 确定颗粒周围 refined 网格的区域
void cfdemCloudMix::setInterface(volScalarField& interface, const double scale) const {
  // reset interface
  interface = dimensionedScalar("zero", interface.dimensions(), 0.0);
  for (int index = 0; index < numberOfParticles(); ++index) {
    if (checkCoarseParticle(index)) {
      Foam::vector particlePos = getPosition(index);
      double radius = getRadius(index);
      forAll(mesh_.C(), cellI) {
        double value = voidFractionM().pointInParticle(mesh_.C()[cellI], particlePos, radius, scale);
        if (value <= 0.0) {
          interface[cellI] = value + 1.0;
        }
      }
    }
  }
}

void cfdemCloudMix::calcVelocityCorrection(volScalarField& p, volVectorField& U, volScalarField& phiIB) const {
  if (validCouplingStep_) {
    // set particle velocity
    Foam::vector parPos = Foam::vector::zero;
    Foam::vector lVel = Foam::vector::zero;
    Foam::vector angVel = Foam::vector::zero;
    Foam::vector rVec = Foam::vector::zero;
    Foam::vector parVel = Foam::vector::zero;
    for (int index = 0; index < numberOfParticles(); ++index) {
      if (checkCoarseParticle(index)) {
        parPos = getPosition(index);         // 颗粒中心
        lVel = getVelocity(index);           // 颗粒线速度
        angVel = getAngularVelocity(index);  // 颗粒角速度
        for (int subCell = 0; subCell < particleOverMeshNumber()[index]; ++subCell) {
          int cellID = cellIDs()[index][subCell];
          if (cellID >= 0) {
            for (int i = 0; i < 3; ++i) {
              rVec[i] = U.mesh().C()[cellID][i] - parPos[i];
            }
            // 计算颗粒速度
            parVel = lVel + (angVel ^ rVec);
            double vf = volumeFractions()[index][subCell];
            U[cellID] = (1.0 - vf) * parVel + vf * U[cellID];
          }
        }
      }
    }
    U.correctBoundaryConditions();
    // phiIB correction equation should be fvm::laplacian(phiIB) == fvc::div(U) according to:
    // https://www.researchgate.net/profile/Stefan_Pirker/publication/264439676_Models_algorithms_and_validation_for_opensource_DEM_and_CFD-DEM/links/56af5af108ae28588c62fd16.pdf
    fvScalarMatrix phiIBEqn(fvm::laplacian(phiIB) == fvc::div(U));
    if (phiIB.needReference()) {
      phiIBEqn.setReference(pRefCell_, pRefValue_);
    }
    phiIBEqn.solve();
    U = U - fvc::grad(phiIB);
    U.correctBoundaryConditions();
    // correct the pressure as well
    p = p + phiIB / U.mesh().time().deltaT();
    p.correctBoundaryConditions();
  }
}

/*!
 * \brief 更新函数
 * \param U          <[in] 流体速度场
 * \param voidF      <[in, out] 小颗粒空隙率场
 * \param volumeF    <[in, out] 大颗粒空隙率场
 * \param Us         <[in, out] 局部平均速度场
 * \param Ksl        <[in, out] 动量交换场
 * \param interface  <[in, out] 界面场
 */
void cfdemCloudMix::evolve(volVectorField& U, volScalarField& voidF, volScalarField& volumeF, volVectorField& Us,
                           volScalarField& Ksl, volScalarField& interface) {
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
    // 根据 DEM data，设置 interface 并更新网格
    updateMesh(interface);
    // 获取颗粒中心所在的网格索引
    locateM().findCell(parCloud_.findCellIDs());
    // 获取颗粒覆盖当前处理器上的某一个网格索引，如果索引为 -1，则表示颗粒不覆盖当前处理器
    // findMpiCellIDs() 作用于所有颗粒，因为其用于计算颗粒尺度
    locateM().findMpiCell(parCloud_.findMpiCellIDs());
    // 计算颗粒尺度
    voidFractionM().getDimensionRatiosForMix(parCloud_.dimensionRatios());
    // 计算扩展颗粒覆盖当前处理器上的某一个网格索引
    // 必须位于计算颗粒尺度之后，因为需要判断颗粒是否为 middle
    // findExpandedCell() 只会作用于 middle 颗粒
    locateM().findExpandedCell(parCloud_.findExpandedCellIDs(), expandedCellScale());
    // 计算颗粒空隙率
    voidFractionM().setVoidFraction();
    voidF = voidFractionM().voidFractionNext();
    volumeF = voidFractionM().volumeFractionNext();
    // 计算 ddtVoidFraction_
    calcDDtVoidFraction(voidF);
    // 计算局部平局颗粒速度场
    averagingM().setVectorFieldAverage(averagingM().UsNext(), averagingM().UsWeightField(), velocities(),
                                       particleWeights());
    Us = averagingM().UsNext();
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
  } else {
    Info << __func__ << ": no need couple, just update mesh..." << endl;
    // 使用上一个耦合时间步的 DEM data 更新 mesh
    // 虽然 DEM data 没有改变，但是这里仍然需要 update mesh，因为 mesh 需要根据自身的 maxRefinement 等属性更新
    updateMesh(interface);
  }
  Info << __func__ << " - done\n" << endl;
}

}  // namespace Foam
