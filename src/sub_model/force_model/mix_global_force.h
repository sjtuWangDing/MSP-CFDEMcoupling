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

#ifndef __MIX_GLOBAL_FORCE_H__
#define __MIX_GLOBAL_FORCE_H__

#include <unordered_set>
#include "./global_force.h"

namespace Foam {

template <typename DType, typename TensorType>
struct fieldRefine {
  template <bool fieldIsNotVF = true>
  static void op(const TensorType& data, const DType& value, double cellV, double GaussCore, double vf);
  //! \brief 子处理器节点算子
  static DType subProcOp();
  //! \brief 主处理器节点算子
  static DType rootProcOp(const std::vector<TensorType> dataVec);
};

template <>
struct fieldRefine<Foam::vector, base::CDTensor1> {
  //! \brief 计算累计数据
  template <bool>
  static void op(const base::CDTensor1& data, const Foam::vector& value, double cellV, double GaussCore, double vf) {
    // 计算累计物理量
    double core = cellV * GaussCore * vf;
    for (int i = 0; i < 3; ++i) {
      data[i] += core * value[i];
    }
    // 计算累计因数
    data[3] += core;
  }
  //! \brief 子处理器节点直接返回
  static Foam::vector subProcOp() { return Foam::vector::zero; }
  //! \brief 主处理器节点计算 refined field
  static Foam::vector rootProcOp(const std::vector<base::CDTensor1>& dataVec) {
    Foam::vector res = Foam::vector::zero;
    double sumCore = 0.0;
    for (const auto& data : dataVec) {
      for (int i = 0; i < 3; ++i) {
        res[i] += data[i];
      }
      sumCore += data[3];
    }
    return res / sumCore;
  }
};

template <>
struct fieldRefine<Foam::scalar, base::CDTensor1> {
  //! \brief 计算累计数据
  template <bool fieldIsNotVF = true>
  static void op(const base::CDTensor1& data, const Foam::scalar& value, double cellV, double GaussCore = 1.0,
                 double vf = 1.0) {
    if (fieldIsNotVF) {
      double core = cellV * GaussCore * vf;
      // 计算累计物理量
      data[0] += value * core;
      // 计算累计因数
      data[1] += core;
    } else {
      // 计算累计物理量
      data[0] += value * cellV;
      // 计算累计因数
      data[1] += cellV;
    }
  }
  //! \brief 子处理器节点直接返回
  static Foam::scalar subProcOp() { return 0.0; }
  //! \brief 主处理器节点计算 refined data
  static Foam::scalar rootProcOp(const std::vector<base::CDTensor1>& dataVec) {
    // 由主节点计算 refined data
    double res = 0.0;
    double sumCore = 0.0;
    for (const auto& data : dataVec) {
      res += data[0];
      sumCore += data[1];
    }
    return res / sumCore;
  }
};

class mixGlobalForce : public globalForce {
 public:
  cfdemTypeName("mixGlobalForce");

  cfdemDefineNewFunctionAdder(globalForce, mixGlobalForce);

  //! \brief Constructor
  mixGlobalForce(cfdemCloud& cloud);

  //! \brief Destructor
  ~mixGlobalForce();

  //! \brief 每一次耦合中，在 set force 前执行
  void initBeforeSetForce();

  //! \brief 每一次耦合中，在 set force 后执行
  void endAfterSetForce();

  //! \brief 获取颗粒处背景流体速度
  Foam::vector getBackgroundUfluid(const int index) const;

  //! \brief 获取颗粒处背景空隙率
  double getBackgroundVoidFraction(const int index) const;

  //! \brief 获取颗粒处背景流体的 ddtU
  Foam::vector getBackgroundDDtU(const int index) const;

 private:
  /*!
   * \brief 计算颗粒 index 处的背景物理场量
   * \tparam fieldIsNotVF 背景物理场量是否是空隙率场，空隙率场的计算与其他场的方法不同
   * \tparam NDim 基本数据的个数，Eg: for Foam::vector, NDim = 3
   * \tparam FieldType 场类型
   * \tparam DType 场元素类型
   * \tparam BDType 场元素基本数据类型
   */
  template <bool fieldIsNotVF = true, int NDim = Foam::vector::dim, typename FieldType = Foam::volVectorField,
            typename DType = Foam::vector, typename BDType = Foam::scalar>
  DType getBackgroundFieldValue(int index, const FieldType& field) const;

 private:
  //! \brief 颗粒覆盖的扩展网格的索引
  //! \brief map: 颗粒索引 --> 扩展网格索引
  std::unordered_map<int, std::unordered_set<int>> expandedCellMap_;

  //! \brief 背景流体速度
  //! \note map: 颗粒索引 --> 背景流体速度
  std::unordered_map<int, Foam::vector> backgroundUfluidMap_;

  //! \brief 背景流体空隙率
  //! \note map: 颗粒索引 --> 背景流体空隙率
  std::unordered_map<int, double> backgroundVoidFractionMap_;

  //! \brief 背景流体的 ddtU
  //! \note map: 颗粒索引 --> 背景流体的ddtU
  std::unordered_map<int, Foam::vector> backgroundDDtUMap_;
};

/*!
 * \brief 计算颗粒 index 处的背景物理场量
 * \tparam fieldIsNotVF 背景物理场量是否是空隙率场，空隙率场的计算与其他场的方法不同
 * \tparam NDim 基本数据的个数，Eg: for Foam::vector, NDim = 3
 * \tparam FieldType 场类型
 * \tparam DType 场元素类型
 * \tparam BDType 场元素基本数据类型
 */
template <bool fieldIsNotVF /*= true */, int NDim /* = Foam::vector::dim*/,
          typename FieldType /* = Foam::volVectorField */, typename DType /* = Foam::vector */,
          typename BDType /* = Foam::scalar*/
          >
DType mixGlobalForce::getBackgroundFieldValue(int index, const FieldType& field) const {
  typedef base::Tensor<1, BDType> TensorType;
  // 数据集合，data[0 ~ NDim - 1] 为累计数据和，data[NDim] 为累计因数
  TensorType data(base::makeShape1(NDim + 1), 0.0);
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
          gcore = GaussCore(particlePos, cellPos, radius, cloud_.expandedCellScale());
          // 计算累计数据
          fieldRefine<DType, TensorType>::template op<fieldIsNotVF>(data, field[cellID], cellV, gcore,
                                                                    voidFraction_[cellID]);
        }
      }
    }
  }
  base::MPI_Barrier();
  DType res = DType();                                 // 计算结果
  int procId = base::procId();                         // 处理器编号
  int numProc = base::numProc();                       // 处理器数量
  int rootProc = cloud_.particleRootProcIDs()[index];  // 主节点编号，即颗粒所在的处理器编号
  int tag = 100;
  // 非主节点
  if (rootProc != procId) {
    // 发送 data 给主节点(非阻塞)
    MPI_Request request;
    base::MPI_Isend(data, rootProc, tag, &request);
    res = fieldRefine<DType, TensorType>::subProcOp();
  }
  // 主节点
  if (rootProc == procId) {
    std::vector<MPI_Request> rVec(numProc - 1);
    std::vector<MPI_Status> sVec(numProc - 1);
    std::vector<TensorType> dataVec;
    // 接收其他节点的 data
    for (int inode = 0; inode < numProc; ++inode) {
      if (inode == rootProc) {
        dataVec.emplace_back(base::makeShape1(NDim + 1), 0.0);
        base::copyTensor(data, dataVec[inode]);
      } else {
        dataVec.emplace_back(base::makeShape1(NDim + 1), 0.0);
        base::MPI_Irecv(dataVec[inode], inode, tag, rVec.data() + (inode < rootProc ? inode : inode - 1));
      }
    }
    // 主节点等待 Irecv 执行完成
    MPI_Waitall(numProc - 1, rVec.data(), sVec.data());
    // 由主节点计算 refined data
    res = fieldRefine<DType, TensorType>::rootProcOp(dataVec);
  }
  // 所有处理器同步后返回结果
  base::MPI_Barrier();
  return res;
}

}  // namespace Foam

#endif  // __MIX_GLOBAL_FORCE_H__
