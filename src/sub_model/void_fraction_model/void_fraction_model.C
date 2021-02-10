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
  demoModel
\*---------------------------------------------------------------------------*/

#include "./void_fraction_model.h"
#include "sub_model/data_exchange_model/data_exchange_model.h"

namespace Foam {

cfdemDefineTypeName(voidFractionModel);

cfdemDefineNewFunctionMap(voidFractionModel);

cfdemDefineConstructNewFunctionMap(voidFractionModel);

cfdemDefineDestroyNewFunctionMap(voidFractionModel);

cfdmeDefineBaseTypeNew(autoPtr, voidFractionModel, (cfdemCloud & cloud, const dictionary& dict), dict, (cloud));

//! \brief Constructor
voidFractionModel::voidFractionModel(cfdemCloud& cloud)
    : cloud_(cloud),
      weight_(1.0),
      porosity_(1.0),
      voidFractionPrev_(IOobject("voidFractionPrev", cloud.mesh().time().timeName(), cloud.mesh(),
                                 IOobject::READ_IF_PRESENT,  // or MUST_READ,
                                 IOobject::AUTO_WRITE),
                        // cloud.mesh().lookupObject<volScalarField>("voidFraction")
                        cloud.mesh(), dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)),
      voidFractionNext_(IOobject("voidFractionNext", cloud.mesh().time().timeName(), cloud.mesh(),
                                 IOobject::READ_IF_PRESENT,  // or MUST_READ,
                                 IOobject::AUTO_WRITE),
                        // cloud.mesh().lookupObject<volScalarField>("voidFraction")
                        cloud.mesh(), dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)),
      volumeFractionPrev_(IOobject("volumeFractionPrev", cloud.mesh().time().timeName(), cloud.mesh(),
                                   IOobject::READ_IF_PRESENT,  // or MUST_READ,
                                   IOobject::AUTO_WRITE),
                          // cloud.mesh().lookupObject<volScalarField>("volumeFraction")
                          cloud.mesh(), dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)),
      volumeFractionNext_(IOobject("volumeFractionNext", cloud.mesh().time().timeName(), cloud.mesh(),
                                   IOobject::READ_IF_PRESENT,  // or MUST_READ,
                                   IOobject::AUTO_WRITE),
                          // cloud.mesh().lookupObject<volScalarField>("volumeFraction")
                          cloud.mesh(), dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)),
      maxCellsNumPerFineParticle_(1),
      maxCellsNumPerMiddleParticle_(1),
      maxCellsNumPerCoarseParticle_(1) {}

//! \brief Destructor
voidFractionModel::~voidFractionModel() {}

//! \brief 计算颗粒尺寸与其周围网格平均尺寸的比值
void voidFractionModel::getDimensionRatios(const base::CDTensor1& dimensionRatios,
                                           const double scale /* = 1.0 */) const {
  base::fillTensor(cloud_.dimensionRatios(), -1.0);
  std::unordered_set<int> set;  // 颗粒扩展网格的集合
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    set.clear();
    int findCellID = cloud_.findCellIDs()[index];
    if (findCellID < 0) {
      continue;
    }
    // 获取颗粒半径
    double radius = cloud_.getRadius(index);
    // 获取颗粒中心位置
    Foam::vector particlePos = cloud_.getPosition(index);
    // 构建集合
    buildExpandedCellSet(set, findCellID, particlePos, radius, 1.0);
    // 计算颗粒周围网格的平均尺寸
    double sumVCell = 0.0;
    for (auto cellID : set) {
      sumVCell += cloud_.mesh().V()[cellID];
    }
    // 计算颗粒尺寸与颗粒中心所在网格尺寸的比值
    double ratio = pow(sumVCell / static_cast<double>(set.size()), 1.0 / 3.0) / (2.0 * radius);
    cloud_.dimensionRatios()[index] = ratio;
  }  // end loop of particles
}

//! \brief 计算颗粒尺寸与其周围网格平均尺寸的比值
void voidFractionModel::getDimensionRatiosForMix(const double scale /* = 1.0 */) const {
  // reset dimensionRatios
  base::fillTensor(cloud_.dimensionRatios(), -1.0);
  int number = cloud_.numberOfParticles();                      // 颗粒总数
  double radius = 0.0;                                          // 颗粒半径
  Foam::vector particlePos = Foam::vector::zero;                // 颗粒中心坐标
  base::CITensor1 sumCellsNumber(base::makeShape1(number), 0);  // 颗粒覆盖当前求解器的网格数量
  base::CDTensor1 sumCellsV(base::makeShape1(number), 0.0);     // 颗粒覆盖当前求解器的网格总体积
  std::unordered_set<int> set;                                  // 颗粒扩展网格的集合
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    set.clear();
    radius = cloud_.getRadius(index);
    particlePos = cloud_.getPosition(index);
    int findMpiCellID = cloud_.findMpiCellIDs()[index];
    if (findMpiCellID < 0) {
      continue;
    }
    // 构建集合
    buildExpandedCellSet(set, findMpiCellID, particlePos, radius, 1.0);
    // 计算颗粒周围网格的平均尺寸
    double sumV = 0.0;
    for (auto cellID : set) {
      sumV += cloud_.mesh().V()[cellID];
    }
    sumCellsNumber[index] = static_cast<int>(set.size());
    sumCellsV[index] = sumV;
  }
  // 主节点汇总其他节点的 dimension ratio 等信息
  int procId = base::procId();
  int numProc = base::numProc();
  int nTag = 100;
  int vTag = 101;
  if (0 != procId) {
    MPI_Request rq1, rq2;
    // 发送 sumCellsNumber 给主节点( MPI_Isend 非阻塞)
    base::MPI_Isend(sumCellsNumber, 0, nTag, &rq1);
    base::MPI_Isend(sumCellsV, 0, vTag, &rq2);
  }
  if (0 == procId) {
    std::vector<MPI_Request> rVec(numProc);
    std::vector<MPI_Status> sVec(numProc);
    std::vector<base::CITensor1> sumCellsNumberVec;
    std::vector<base::CDTensor1> sumCellsVVec;
    // 接收其他节点的 sumCellsNumber 信息
    for (int inode = 1; inode < numProc; ++inode) {
      sumCellsNumberVec.emplace_back(base::makeShape1(number), 0);
      base::MPI_Irecv(sumCellsNumberVec[inode - 1], inode, nTag, rVec.data() + inode);
    }
    // 主节点等待 Irecv 执行完成
    MPI_Waitall(numProc - 1, rVec.data() + 1, sVec.data() + 1);
    // 接收其他节点的 sumCellsV 信息
    for (int inode = 1; inode < numProc; ++inode) {
      sumCellsVVec.emplace_back(base::makeShape1(number), 0.0);
      base::MPI_Irecv(sumCellsVVec[inode - 1], inode, vTag, rVec.data() + inode);
    }
    // 主节点等待 Irecv 执行完成
    MPI_Waitall(numProc - 1, rVec.data() + 1, sVec.data() + 1);
    // 由主节点计算 dimension ratio
    for (int index = 0; index < number; ++index) {
      radius = cloud_.getRadius(index);
      int Nmesh = 0;
      double Vmesh = 0;
      for (int j = 0; j < numProc; ++j) {
        if (0 == j) {
          Nmesh += std::max(sumCellsNumber[index], 0);
          Vmesh += std::max(sumCellsV[index], 0.0);
        } else {
          Nmesh += std::max(sumCellsNumberVec[j - 1][index], 0);
          Vmesh += std::max(sumCellsVVec[j - 1][index], 0.0);
        }
      }
      CHECK_GT(Nmesh, 0) << __func__ << ": particle " << index
                         << " cover no mpi cell with position = " << cloud_.positions()[index]
                         << " and radius = " << cloud_.getRadius(index);
      double ratio = pow(Vmesh / Nmesh, 1.0 / 3.0) / (2.0 * radius);
      cloud_.dimensionRatios()[index] = ratio;
    }
  }
  // 主节点广播 dimensionRatios
  base::MPI_Bcast(cloud_.dimensionRatios(), 0);
  base::MPI_Barrier();
  base::MPI_Info("voidFractionModel: getDimensionRatiosForMix - done", true);
}

/*!
 * \brief 构建颗粒覆盖的所有网格的集合
 * \param set         <[in, out] 需要构建的集合
 * \param cellID      <[in] 递归循环中要检索网格编号
 * \param particlePos <[in] 颗粒中心位置
 * \param radius      <[in] 颗粒半径
 * \param scale       <[in] 颗粒半径扩大系数
 */
void voidFractionModel::buildExpandedCellSet(std::unordered_set<int>& set, const int cellID,
                                             const Foam::vector& particlePos, const double radius,
                                             const double scale) const {
  // 如果搜索到网格的边界，则结束搜索
  if (cellID < 0) {
    return;
  }
  // 将 cellID 插入到哈希集合中
  set.insert(cellID);
  // 获取 cellID 网格的所有 neighbour cell 的链表
  const labelList& nc = cloud_.mesh().cellCells()[cellID];
  // 网格中心坐标
  Foam::vector neighbourPos = Foam::vector::zero;
  // 遍历链表
  for (int i = 0; i < nc.size(); ++i) {
    int neighbourCellID = nc[i];
    neighbourPos = cloud_.mesh().C()[neighbourCellID];
    // 在集合中没有 neighbour 网格的索引
    if (set.end() == set.find(neighbourCellID)) {
      if (pointInParticle(neighbourPos, particlePos, radius, scale) < 0.0) {
        // neighbour 网格中心在颗粒中，以 neighbour 为中心递归构建集合
        buildExpandedCellSet(set, neighbourCellID, particlePos, radius, scale);
      } else {
        bool vertexPointInParticle = false;
        // 遍历网格的角点
        const labelList& vertexPoints = cloud_.mesh().cellPoints()[neighbourCellID];
        // 遍历当前网格的所有角点
        forAll(vertexPoints, i) {
          // 获取第 i 角点坐标
          vector vertexPosition = cloud_.mesh().points()[vertexPoints[i]];
          // 判断角点是否在颗粒中
          double fv = pointInParticle(vertexPosition, particlePos, radius, scale);
          if (fv < 0.0) {
            vertexPointInParticle = true;
            break;
          }
        }
        if (vertexPointInParticle) {
          // 如果有一个角点在颗粒中，则以 neighbour 为中心递归构建集合
          buildExpandedCellSet(set, neighbourCellID, particlePos, radius, scale);
        }
      }
    }
  }  // end loop of neighbour list
}

}  // namespace Foam
