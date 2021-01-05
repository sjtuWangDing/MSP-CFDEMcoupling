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
  cfdemCloudIBOpti derived from cfdemCloud

Class
  Foam::cfdemCloudIBOpti
\*---------------------------------------------------------------------------*/

#include <mutex>
#include "cloud/cfdem_cloud_IB_opti.h"
#include "dynamicFvMesh.H"
#include "sub_model/data_exchange_model/data_exchange_model.h"
#include "sub_model/force_model/force_model.h"

namespace Foam {

//! \brief Constructed from mesh
cfdemCloudIBOpti::cfdemCloudIBOpti(const fvMesh& mesh)
    : cfdemCloud(mesh),
      pRefCell_(readLabel(mesh.solutionDict().subDict("PISO").lookup("pRefCell"))),
      pRefValue_(readScalar(mesh.solutionDict().subDict("PISO").lookup("pRefValue"))) {}

//! \brief Destructor
cfdemCloudIBOpti::~cfdemCloudIBOpti() {
  dataExchangeM().free(parCloud_.radii(), parCloud_.radiiPtr());
  dataExchangeM().free(parCloud_.positions(), parCloud_.positionsPtr());
  dataExchangeM().free(parCloud_.velocities(), parCloud_.velocitiesPtr());
  dataExchangeM().free(parCloud_.angularVelocities(), parCloud_.angularVelocitiesPtr());
}

//! \brief 重新分配内存
void cfdemCloudIBOpti::reallocate() {
  int number = numberOfParticles();
  // allocate memory of data exchanged with liggghts
  dataExchangeM().realloc(parCloud_.radii(), base::makeShape1(number), parCloud_.radiiPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.positions(), base::makeShape2(number, 3), parCloud_.positionsPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.velocities(), base::makeShape2(number, 3), parCloud_.velocitiesPtr(), 0.0);
  dataExchangeM().realloc(parCloud_.angularVelocities(), base::makeShape2(number, 3), parCloud_.angularVelocitiesPtr(),
                          0.0);
  dataExchangeM().realloc(parCloud_.DEMForces(), base::makeShape2(number, 3), parCloud_.DEMForcesPtr(), 0.0);
  // allocate memory of data not exchanged with liggghts
  parCloud_.particleOverMeshNumber() = std::move(base::CITensor1(base::makeShape1(number), 0));
  parCloud_.findCellIDs() = std::move(base::CITensor1(base::makeShape1(number), -1));
  parCloud_.dimensionRatios() = std::move(base::CDTensor1(base::makeShape1(number), -1.0));
}

void cfdemCloudIBOpti::getDEMData() {
  dataExchangeM().getData("radius", "scalar-atom", parCloud_.radiiPtr());
  dataExchangeM().getData("x", "vector-atom", parCloud_.positionsPtr());
  dataExchangeM().getData("v", "vector-atom", parCloud_.velocitiesPtr());
  dataExchangeM().getData("omega", "vector-atom", parCloud_.angularVelocitiesPtr());
}

void cfdemCloudIBOpti::giveDEMData() const {
  dataExchangeM().giveData("dragforce", "vector-atom", DEMForcesPtr());
}

/*!
 * \brief 更新网格，如果 mesh 是 Foam::dynamicRefineFvMesh 类型，则更新网格，如果是 Foam::staticFvMesh
 * 或者其他类型，则不更新
 */
void cfdemCloudIBOpti::updateMesh(volScalarField& interface) {
  Info << "update mesh..." << endl;
  fvMesh& mesh = const_cast<fvMesh&>(mesh_);
  try {
    // 采用动态网格
    dynamicFvMesh& dyMesh = dynamic_cast<dynamicFvMesh&>(mesh);
    // 设置 interface
    setInterface(interface);
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
  Info << "update mesh - done" << endl;
  base::MPI_Barrier();
  return;
}

//! \brief 确定颗粒周围 refined 网格的区域
void cfdemCloudIBOpti::setInterface(volScalarField& interface,
                                    const double scale /* = cfdemCloudIBOpti::particleMeshScale_ */) const {
  interface = dimensionedScalar("zero", interface.dimensions(), 0.0);
  for (int index = 0; index < numberOfParticles(); ++index) {
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

void cfdemCloudIBOpti::calcVelocityCorrection(volScalarField& p, volVectorField& U, volScalarField& phiIB,
                                              volScalarField& volumeFraction) {
  if (validCouplingStep_) {
    // set particle velocity
    Foam::vector parPos = Foam::vector::zero;
    Foam::vector lVel = Foam::vector::zero;
    Foam::vector angVel = Foam::vector::zero;
    Foam::vector rVec = Foam::vector::zero;
    Foam::vector parVel = Foam::vector::zero;
    for (int index = 0; index < numberOfParticles(); ++index) {
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
    U.correctBoundaryConditions();
    fvScalarMatrix phiIBEqn(fvm::laplacian(phiIB) == fvc::div(U) + fvc::ddt(volumeFraction));
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
 * \note used for cfdemSolverIB
 * \param volumeFraction  <[in, out] 大颗粒体积分数
 * \param interface       <[in, out] 界面场，用于 dynamic mesh
 */
void cfdemCloudIBOpti::evolve(volScalarField& volumeFraction, volScalarField& interface) {
  Info << "/ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /" << endl;
  // 检查当前流体时间步是否同时也是耦合时间步
  validCouplingStep_ = dataExchangeM().checkValidCouplingStep();
  Info << "time step fraction: " << dataExchangeM().timeStepFraction() << endl;
  if (validCouplingStep_) {
    Info << __func__ << ": coupling..." << endl;
    // 创建用于记录 coupling time step counter
    auto pCounter = std::make_shared<dataExchangeModel::CouplingStepCounter>(dataExchangeM());
    // couple(): run liggghts command and get number of particle
    setNumberOfParticles(dataExchangeM().couple());
    // realloc memory
    reallocate();
    // 获取 DEM data
    getDEMData();
    // 根据 DEM data，设置 interface 并更新网格
    updateMesh(interface);
    // 获取到在当前 processor 上颗粒覆盖的某一个网格编号，如果获取到的网格编号为 -1，则表示颗粒不覆盖当前 processor
    locateM().findCell(parCloud_.findCellIDs());
    // 计算颗粒空隙率
    voidFractionM().setVoidFraction();
    volumeFraction == voidFractionM().volumeFractionNext();
    volumeFraction.correctBoundaryConditions();
    // reset particles forces
    std::fill_n(DEMForces().ptr(), DEMForces().mSize(), 0.0);
    // set particles forces
    for (const auto& ptr : forceModels_) {
      ptr->setForce();
    }
    // write DEM data
    giveDEMData();
    // 输出信息
    printParticleInfo();
  } else {
    Info << __func__ << ": no need couple, just update mesh..." << endl;
    // 使用上一个耦合时间步的 DEM data 更新 mesh
    // 虽然 DEM data 没有改变，但是这里仍然需要 update mesh，因为 mesh 需要根据自身的 maxRefinement 等属性更新
    updateMesh(interface);
  }
  Info << __func__ << " - done\n" << endl;
}

// //! \brief 确定颗粒周围 refined 网格的区域(每个方向的尺寸都是颗粒尺寸的两倍)
// void cfdemCloudIB::setInterface(volScalarField& interface, volScalarField& refineMeshKeepStep) const {
//   // 确保 call_once 函数中的 lambda 表达式只会被执行一次
//   static std::once_flag flag;
//   std::call_once(flag, [&interface, &refineMeshKeepStep]() {
//     // reset interface and refineMeshKeepStep field
//     // only set at first step ?????????????????????????????????????
//     interface = dimensionedScalar("zero", interface.dimensions(), 0.0);
//     refineMeshKeepStep = dimensionedScalar("zero", refineMeshKeepStep.dimensions(), 0.0);
//   });
//   forAll(mesh_.C(), cellI) {
//     // 当前 cellI 是否位于任意一个颗粒中
//     bool cellInParticle = false;
//     bool cellFirstEntryRefineMeshKeepStep = false;
//     // 网格中心
//     Foam::vector cellPos = mesh_.C()[cellI];
//     for (int index = 0; index < numberOfParticles(); ++index) {
//       if (checkPeriodicCells()) {
//         FatalError << "Error: not support periodic check!" << abort(FatalError);
//       }
//       // 判断网格中心是否在 index 颗粒中
//       double value = voidFractionM().pointInParticle(getPosition(index), cellPos, getRadius(index),
//       refineMeshSkin());
//       if (0 == refineMeshKeepInterval()) {
//         if (value <= 0.0) {
//           interface[cellI] = std::max(interface[cellI], value + 1);
//         }
//       } else {  // valid refineMeshKeepInterval
//         if (value <= 0.0) {
//           // cellI 位于 index 颗粒内部
//           cellInParticle = true;
//           // 如果 cellI 网格位于任何一个颗粒内部，则重新设置 refineMeshKeepStep
//           refineMeshKeepStep[cellI] = refineMeshKeepInterval();
//           // 设置 interface
//           interface[cellI] = std::max(interface[cellI], value + 1);
//         } else {
//           if (false == cellInParticle) {
//             // cellI 目前不在任何颗粒中
//             if (refineMeshKeepStep[cellI] > Foam::SMALL && false == cellFirstEntryRefineMeshKeepStep) {
//               // refineMeshKeepStep[cellI] > 0.0，则保持 interFace 值
//               cellFirstEntryRefineMeshKeepStep = true;  // 确保对每一个 cellI 只执行一次
//               refineMeshKeepStep[cellI] -= 1.0;
//             } else {
//               // 设置 interFace 为 0.0
//               interface[cellI] = 0.0;
//             }
//           }
//         }
//       }
//     }  // end of all particles
//   }    // end of all cell in current processor
// }

}  // namespace Foam
