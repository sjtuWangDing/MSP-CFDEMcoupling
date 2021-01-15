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
  cfdemCloudIB derived from cfdemCloud

Class
  Foam::cfdemCloudIB
\*---------------------------------------------------------------------------*/

#include <mutex>
#include "cloud/cfdem_cloud_IB.h"
#include "mpi.h"
#include "sub_model/data_exchange_model/data_exchange_model.h"
#include "sub_model/force_model/force_model.h"
#include "sub_model/liggghts_command_model/liggghts_command_model.h"
#include "sub_model/locate_model/locate_model.h"
#include "sub_model/void_fraction_model/void_fraction_model.h"

namespace Foam {

//! \brief Constructed from mesh
cfdemCloudIB::cfdemCloudIB(const fvMesh& mesh)
    : cfdemCloud(mesh),
      pRefCell_(readLabel(mesh.solutionDict().subDict("PISO").lookup("pRefCell"))),
      pRefValue_(readScalar(mesh.solutionDict().subDict("PISO").lookup("pRefValue"))) {}

//! \brief Destructor
cfdemCloudIB::~cfdemCloudIB() {
  dataExchangeM().free(parCloud_.radii(), parCloud_.radiiPtr());
  dataExchangeM().free(parCloud_.positions(), parCloud_.positionsPtr());
  dataExchangeM().free(parCloud_.velocities(), parCloud_.velocitiesPtr());
  dataExchangeM().free(parCloud_.angularVelocities(), parCloud_.angularVelocitiesPtr());
}

//! \brief 重新分配内存
void cfdemCloudIB::reallocate() {
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

void cfdemCloudIB::getDEMData() {
  dataExchangeM().getData("radius", "scalar-atom", parCloud_.radiiPtr());
  dataExchangeM().getData("x", "vector-atom", parCloud_.positionsPtr());
  dataExchangeM().getData("v", "vector-atom", parCloud_.velocitiesPtr());
  dataExchangeM().getData("omega", "vector-atom", parCloud_.angularVelocitiesPtr());
}

void cfdemCloudIB::giveDEMData() const {
  dataExchangeM().giveData("dragforce", "vector-atom", DEMForcesPtr());
}

void cfdemCloudIB::printParticleInfo() const {
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
    Pout << "  DEMForce[" << index << "]: " << DEMForces()[index][0] << ", " << DEMForces()[index][1] << ", "
         << DEMForces()[index][2] << endl;
  }
  base::MPI_Barrier();
  voidFractionM().printVoidFractionInfo();
}

//! \brief 确定颗粒周围 refined 网格的区域
void cfdemCloudIB::setInterface(volScalarField& interface,
                                const double scale /* = cfdemCloudIB::particleMeshScale_ */) const {
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

void cfdemCloudIB::calcVelocityCorrection(volScalarField& p, volVectorField& U, volScalarField& phiIB,
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
void cfdemCloudIB::evolve(volScalarField& volumeFraction, volScalarField& interface) {
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
    // 获取到在当前 processor 上颗粒覆盖的某一个网格编号，如果获取到的网格编号为 -1，则表示颗粒不覆盖当前 processor
    locateM().findCell(parCloud_.findCellIDs());
    setInterface(interface);
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
  }
  Info << __func__ << " - done\n" << endl;
}

}  // namespace Foam
