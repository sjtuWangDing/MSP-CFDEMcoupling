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

#include "./mix_Basset_force.h"
#include "./mix_Mei_lift_force.h"
#include "./mix_global_force.h"
#include "./mix_virtual_mass_force.h"
#include "cfdem_tools/cfdem_tools.h"
#include "sub_model/void_fraction_model/void_fraction_model.h"

namespace Foam {

cfdemDefineTypeName(mixGlobalForce);

cfdemCreateNewFunctionAdder(globalForce, mixGlobalForce);

//! \brief Constructor
mixGlobalForce::mixGlobalForce(cfdemCloud& cloud) : globalForce(cloud) {}

//! \brief Destructor
mixGlobalForce::~mixGlobalForce() {}

#if 0

//! \brief 每一次耦合中，在 set force 前执行
void mixGlobalForce::initBeforeSetForce() {
  // init data - Important !!!
  // nerver forget !!!
  expandedCellMap_.clear();
  backgroundUfluidMap_.clear();
  backgroundVoidFractionMap_.clear();
  backgroundDDtUMap_.clear();
  backgroundVorticityMap_.clear();

  // 如果使用了 BassetForce or virtualMassForce，则需要计算 ddtU
  bool isUsedVirtualMassForce = cfdemTools::isUsedForceModel(cloud_, mixVirtualMassForce::cTypeName());
  bool isUsedBassetForce = cfdemTools::isUsedForceModel(cloud_, mixBassetForce::cTypeName());
  if (isUsedVirtualMassForce || isUsedBassetForce) {
    // 计算 ddtU field
    ddtU_ = fvc::ddt(U_) + fvc::div(phi_, U_);
  }

  // 如果使用了升力，则需要计算 vorticityField
  bool isUsedMixMeiLiftForce = cfdemTools::isUsedForceModel(cloud_, mixMeiLiftForce::cTypeName());
  if (isUsedMixMeiLiftForce) {
    // 计算 vorticityField
    vorticityField_ = fvc::curl(U_);
  }

  // reset data
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkMiddleParticle(index)) {
      double radius = cloud_.getRadius(index);                       // 颗粒半径
      Foam::vector particlePos = cloud_.getPosition(index);          // 颗粒中心坐标
      int findExpandedCellID = cloud_.findExpandedCellIDs()[index];  // 扩展网格ID
      // 计算颗粒覆盖的扩展网格集合
      expandedCellMap_.insert(std::make_pair(index, std::unordered_set<int>()));
      if (findExpandedCellID >= 0) {
        cloud_.voidFractionM().buildExpandedCellSet(expandedCellMap_[index], findExpandedCellID, particlePos, radius,
                                                    cloud_.expandedCellScale());
      }
      // 计算背景流体速度
      backgroundUfluidMap_.insert(std::make_pair(index, getBackgroundFieldValue(index, U_)));
      // 计算背景流体空隙率
      backgroundVoidFractionMap_.insert(
          std::make_pair(index, getBackgroundFieldValue<false, 1, volScalarField, scalar>(index, voidFraction_)));
      // 计算背景流体的 ddtU
      if (isUsedVirtualMassForce || isUsedBassetForce) {
        backgroundDDtUMap_.insert(std::make_pair(index, getBackgroundFieldValue(index, ddtU_)));
      }
      // 计算背景流体的涡量
      if (isUsedMixMeiLiftForce) {
        backgroundVorticityMap_.insert(std::make_pair(index, getBackgroundFieldValue(index, vorticityField_)));
      }
    }
  }
  base::MPI_Info("mixGlobalForce: initBeforeSetForce - done", verbose_);
}

#else

//! \brief 每一次耦合中，在 set force 前执行
void mixGlobalForce::initBeforeSetForce() {
  // (1) init data - Important !!!
  // nerver forget !!!
  expandedCellMap_.clear();
  backgroundUfluidMap_.clear();
  backgroundVoidFractionMap_.clear();
  backgroundDDtUMap_.clear();
  backgroundVorticityMap_.clear();

  // (2) 如果使用了 BassetForce or virtualMassForce，则需要计算 ddtU
  bool isUsedVirtualMassForce = cfdemTools::isUsedForceModel(cloud_, mixVirtualMassForce::cTypeName());
  bool isUsedBassetForce = cfdemTools::isUsedForceModel(cloud_, mixBassetForce::cTypeName());
  if (isUsedVirtualMassForce || isUsedBassetForce) {
    // 计算 ddtU field
    ddtU_ = fvc::ddt(U_) + fvc::div(phi_, U_);
  }

  // (3) 如果使用了升力，则需要计算 vorticityField
  bool isUsedMixMeiLiftForce = cfdemTools::isUsedForceModel(cloud_, mixMeiLiftForce::cTypeName());
  if (isUsedMixMeiLiftForce) {
    // 计算 vorticityField
    vorticityField_ = fvc::curl(U_);
  }

  // (4) build expanded cell set
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkMiddleParticle(index)) {
      double radius = cloud_.getRadius(index);                       // 颗粒半径
      Foam::vector particlePos = cloud_.getPosition(index);          // 颗粒中心坐标
      int findExpandedCellID = cloud_.findExpandedCellIDs()[index];  // 扩展网格ID
      // 计算颗粒覆盖的扩展网格集合
      expandedCellMap_.insert(std::make_pair(index, std::unordered_set<int>()));
      if (findExpandedCellID >= 0) {
        cloud_.voidFractionM().buildExpandedCellSet(expandedCellMap_[index], findExpandedCellID, particlePos, radius,
                                                    cloud_.expandedCellScale());
      }
    }
  }
  base::MPI_Barrier();

  // (5) 计算背景流体速度
  setBackgroundFieldValue<true, 3, volVectorField, Foam::vector>(U_, backgroundUfluidMap_);
  // (6) 计算背景流体空隙率
  setBackgroundFieldValue<false, 1, volScalarField, Foam::scalar>(voidFraction_, backgroundVoidFractionMap_);
  // (7) 计算背景流体的 ddtU
  if (isUsedVirtualMassForce || isUsedBassetForce) {
    setBackgroundFieldValue<true, 3, volVectorField, Foam::vector>(ddtU_, backgroundDDtUMap_);
  }
  // (8) 计算背景流体的涡量
  if (isUsedMixMeiLiftForce) {
    setBackgroundFieldValue<true, 3, volVectorField, Foam::vector>(vorticityField_, backgroundVorticityMap_);
  }
  base::MPI_Info("mixGlobalForce: initBeforeSetForce - done", verbose_);

#if 0

  int number = cloud_.numberOfParticles();
  base::CDTensor2 dataTensor(base::makeShape2(number, 4), 0.0);
  for (int index = 0; index < number; ++index) {
    int findExpandedCellID = cloud_.findExpandedCellIDs()[index];
    if (findExpandedCellID >= 0) {
      Foam::vector particlePos = cloud_.getPosition(index);  // 颗粒中心坐标
      Foam::vector cellPos = Foam::vector::zero;             // 网格中心坐标
      double radius = cloud_.getRadius(index);               // 颗粒半径
      double gcore = 0.0;                                    // 高斯核
      double cellV = 0.0;                                    // 网格体积
      auto iter = expandedCellMap_.find(index);
      if (expandedCellMap_.end() != iter) {
        const std::unordered_set<int>& set = iter->second;  // 颗粒扩展网格的集合
        for (int cellID : set) {
          if (cellID >= 0) {  // cell found
            cellPos = cloud_.mesh().C()[cellID];
            cellV = cloud_.mesh().V()[cellID];
            // 计算高斯核
            gcore = GaussCore(particlePos, cellPos, radius);
            // 计算累计数据
            fieldRefine<Foam::vector, base::CDTensor1>::op<true>(dataTensor[index], U_[cellID], cellV, gcore,
                                                                 voidFraction_[cellID]);
          }
        }
      }
    }
  }
  base::MPI_Barrier();
  int procId = base::procId();    // 处理器编号
  int numProc = base::numProc();  // 处理器数量
  int tag = 100;
  if (0 != procId) {
    MPI_Request request;
    MPI_Status status;
    // 发送 dataTensor 给主节点(非阻塞)
    base::MPI_Isend(dataTensor, 0, tag, &request);
    MPI_Wait(&request, &status);
  }
  if (0 == procId) {
    std::vector<MPI_Request> rVec(numProc);
    std::vector<MPI_Status> sVec(numProc);
    std::vector<base::CDTensor2> dataTensorVec;
    // 接收其他节点的 dataTensor
    for (int inode = 0; inode < numProc; ++inode) {
      if (inode == 0) {
        dataTensorVec.emplace_back(base::makeShape2(number, 4), 0.0);
        base::copyTensor(dataTensor, dataTensorVec[inode]);
      } else {
        dataTensorVec.emplace_back(base::makeShape2(number, 4), 0.0);
        base::MPI_Irecv(dataTensorVec[inode], inode, tag, rVec.data() + inode);
      }
    }
    // 主节点等待 Irecv 执行完成
    MPI_Waitall(numProc - 1, rVec.data() + 1, sVec.data() + 1);
    // 主节点重置 dataTensor
    std::memset(dataTensor.ptr(), dataTensor.mSize(), 0.0);
    for (int index = 0; index < number; ++index) {
      for (int i = 0; i < 4; ++i) {
        dataTensor[index][i] = 0.0;
      }
      if (index == 0) {
        Info << dataTensor[0][0] << ", " << dataTensor[0][1] << ", " << dataTensor[0][2] << ", " << dataTensor[0][3]
             << ", " << endl;
      }
      for (int inode = 0; inode < numProc; ++inode) {
        for (int i = 0; i < 3; ++i) {
          dataTensor[index][i] += dataTensorVec[inode][index][i];
        }
        dataTensor[index][3] += dataTensorVec[inode][index][3];
      }
      for (int i = 0; i < 3; ++i) {
        dataTensor[index][i] /= dataTensor[index][3];
      }
      if (index == 0) {
        Info << dataTensor[0][0] << ", " << dataTensor[0][1] << ", " << dataTensor[0][2] << ", " << dataTensor[0][3]
             << ", " << endl;
      }
    }
  }

  // 主节点广播 dataTensor
  MPI_Bcast(dataTensor.ptr(), dataTensor.mSize(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  base::MPI_Barrier();

  for (int index = 0; index < number; ++index) {
    // 计算背景流体速度
    int rootProc = cloud_.particleRootProcIDs()[index];  // 主节点编号，即颗粒所在的处理器编号
    if (rootProc == procId) {
      backgroundUfluidMap_.insert(
          std::make_pair(index, Foam::vector(dataTensor[index][0], dataTensor[index][1], dataTensor[index][2])));
    } else {
      backgroundUfluidMap_.insert(std::make_pair(index, Foam::vector::zero));
    }
  }
#endif
}

#endif

//! \brief 每一次耦合中，在 set force 后执行
void mixGlobalForce::endAfterSetForce() {
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // (1) 更新颗粒速度，所有处理器都必须更新，当颗粒从一个计算域运动到另一个计算域，新的计算域必须有颗粒的 prevU
    updatePrevParticleVelMap(index);
    // (2) 更新 ddtUrHistoryMap_，所有处理器都必须更新
    if (cloud_.checkMiddleParticle(index) && cfdemTools::isUsedForceModel(cloud_, mixBassetForce::cTypeName())) {
      updateDDtUrHistoryMap(index);
    }
  }
  base::MPI_Info("mixGlobalForce: endAfterSetForce - done", verbose_);
}

//! \brief 获取颗粒处背景流体速度
Foam::vector mixGlobalForce::getBackgroundUfluid(const int index) const {
  CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
  auto iter = backgroundUfluidMap_.find(index);
  if (backgroundUfluidMap_.end() != iter) {
    return iter->second;
  }
  return Foam::vector::zero;
}

//! \brief 获取颗粒处背景空隙率
double mixGlobalForce::getBackgroundVoidFraction(const int index) const {
  CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
  auto iter = backgroundVoidFractionMap_.find(index);
  if (backgroundVoidFractionMap_.end() != iter) {
    return iter->second;
  }
  return 1.0;
}

//! \brief 获取颗粒处背景流体的 ddtU
Foam::vector mixGlobalForce::getBackgroundDDtU(const int index) const {
  CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
  auto iter = backgroundDDtUMap_.find(index);
  if (backgroundDDtUMap_.end() != iter) {
    return iter->second;
  }
  return Foam::vector::zero;
}

//! \brief 获取颗粒处背景流体的涡量
Foam::vector mixGlobalForce::getBackgroundVorticity(const int index) const {
  CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
  auto iter = backgroundVorticityMap_.find(index);
  if (backgroundVorticityMap_.end() != iter) {
    return iter->second;
  }
  return Foam::vector::zero;
}

}  // namespace Foam