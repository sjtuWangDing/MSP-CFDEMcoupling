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
  The locateModel “engine” locates the CFD cell and cellID corresponding
  to a given position.

Syntax
  locateModel engine;
  engineProps
  {
    treeSearch switch1;
  }

Class
  engineSearch
\*---------------------------------------------------------------------------*/

#include "./engine_search.h"

namespace Foam {

cfdemDefineTypeName(engineSearch);

cfdemCreateNewFunctionAdder(locateModel, engineSearch);

//! \brief Constructor
engineSearch::engineSearch(cfdemCloud& cloud, const std::string& derivedTypeName)
    : locateModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(derivedTypeName + "Props")),
      treeSearch_(subPropsDict_.lookupOrDefault<bool>("treeSearch", true)),
#if defined(version30)
      searchEngine_(cloud.mesh(), polyMesh::FACE_PLANES)
#elif defined(version21)
      searchEngine_(cloud.mesh(), polyMesh::FACEPLANES)
#elif defined(version16ext)
      searchEngine_(cloud.mesh(), false)
#endif
{
}

//! \brief Destructor
engineSearch::~engineSearch() {}

/*!
 * \brief use search engine to get cell id of particle center
 * \param findCellIDs 颗粒覆盖网格的编号
 */
void engineSearch::findCell(const base::CITensor1& findCellIDs) const {
  // 初始化 findCellIDs
  std::fill_n(findCellIDs.ptr(), findCellIDs.mSize(), -1);
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // arg1: location vector
    // arg2: old cell id
    // arg3: whether use octree search
    findCellIDs[index] = searchEngine_.findCell(cloud_.getPosition(index), -1, treeSearch_);
  }
  uniqueFindCell(findCellIDs);
  // 检查颗粒中心所在的处理器上，findCellID 是否 >= 0
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    int rootProc = cloud_.particleRootProcIDs()[index];
    if (rootProc == base::procId()) {
      Pout << __func__ << ": particle " << index << " is on proc " << rootProc << endl;
      CHECK_GE(findCellIDs[index], 0) << ": findCellIDs[" << index << "] not GE -1";
    } else {
      CHECK_EQ(findCellIDs[index], -1) << ": findCellIDs[" << index << "] not EQ -1";
    }
    base::MPI_Barrier();
  }
}

//! \brief 每个颗粒中心只可能位于一个网格中，但是如果颗粒位于处理器边界上，则颗粒会被边界两边的处理器都捕获到
//!   所以该函数会确保所有颗粒只被一个处理器捕获到
void engineSearch::uniqueFindCell(const base::CITensor1& findCellIDs) const {
  base::MPI_Barrier();
  // 主节点汇总其他节点的 findCellIDs
  int number = cloud_.numberOfParticles();
  int procId = base::procId();
  int numProc = base::numProc();
  int tag1 = 100;
  int tag2 = 101;
  if (0 != procId) {
    MPI_Request request1, request2, request3;
    MPI_Status status1, status2;
    // 发送 findCellIDs 给主节点(非阻塞)
    base::MPI_Isend(findCellIDs, 0, tag1, &request1);
    // 从主节点接收 findCellIDs
    base::MPI_Irecv(findCellIDs, 0, tag1, &request2);
    // 从主节点接收 particleRootProcIDs
    base::MPI_Irecv(cloud_.particleRootProcIDs(), 0, tag2, &request3);
    MPI_Wait(&request2, &status1);
    MPI_Wait(&request3, &status2);
  }
  if (0 == procId) {
    std::vector<MPI_Request> rVec(numProc);
    std::vector<MPI_Status> sVec(numProc);
    std::vector<base::CITensor1> findCellIDsVec;
    // 将主节点的数据拷贝到 findCellIDsVec 中
    findCellIDsVec.emplace_back(base::makeShape1(number), -1);
    base::copyTensor(findCellIDs, findCellIDsVec[0]);
    // 接收其他节点的 findCellIDs 信息
    for (int inode = 1; inode < numProc; ++inode) {
      findCellIDsVec.emplace_back(base::makeShape1(number), 0);
      base::MPI_Irecv(findCellIDsVec[inode], inode, tag1, rVec.data() + inode);
    }
    // 主节点等待 Irecv 执行完成
    MPI_Waitall(numProc - 1, rVec.data() + 1, sVec.data() + 1);
    // 由主节点统计 findCellIDs
    for (int index = 0; index < number; ++index) {
      for (int inode = 0; inode < numProc; ++inode) {
        if (findCellIDsVec[inode][index] == -1) {
          continue;
        } else if (findCellIDsVec[inode][index] >= 0) {
          // index 颗粒中心位于编号为 inode 的处理器上
          cloud_.particleRootProcIDs()[index] = inode;
          for (int j = inode + 1; j < numProc; ++j) {
            findCellIDsVec[j][index] = -1;
          }
          break;
        }
      }
    }
    // 由主节点将 findCellIDsVec 传递给子节点
    for (int inode = 1; inode < numProc; ++inode) {
      base::MPI_Isend(findCellIDsVec[inode], inode, tag1, rVec.data() + inode);
      base::MPI_Isend(cloud_.particleRootProcIDs(), inode, tag2, rVec.data() + inode);
    }
  }
  base::MPI_Barrier();
}

/*!
 * \brief use search engine to get cell id of certain vector
 * \param position 颗粒中心的位置
 * \param oldCellID old cell ID
 */
label engineSearch::findSingleCell(const Foam::vector& position, label oldCellID) const {
  return searchEngine_.findCell(position, oldCellID, treeSearch_);
}

}  // namespace Foam
