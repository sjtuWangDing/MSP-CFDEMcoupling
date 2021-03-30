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
  static void op(const TensorType& data, const DType& value, double cellV, double GaussCore, double voidF,
                 double volumeF);
  //! \brief 子处理器节点算子
  static DType subProcOp();
  //! \brief 主处理器节点算子
  static DType rootProcOp(const std::vector<TensorType> dataVec);
  static DType rootProcOp(const TensorType& dataT);
};

template <>
struct fieldRefine<Foam::vector, base::CDTensor1> {
  //! \brief 计算累计数据
  template <bool>
  static void op(const base::CDTensor1& data, const Foam::vector& value, double cellV, double GaussCore, double voidF,
                 double volumeF) {
    // 计算累计物理量
    double core = cellV * GaussCore * voidF * volumeF;
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
  static Foam::vector rootProcOp(const base::CDTensor1& dataT) { return Foam::vector(dataT[0], dataT[1], dataT[2]); }
};

template <>
struct fieldRefine<Foam::scalar, base::CDTensor1> {
  //! \brief 计算累计数据
  template <bool fieldIsNotVF = true>
  static void op(const base::CDTensor1& data, const Foam::scalar& value, double cellV, double GaussCore, double voidF,
                 double volumeF) {
    if (fieldIsNotVF) {
      double core = cellV * GaussCore * voidF * volumeF;
      // 计算累计物理量
      data[0] += value * core;
      // 计算累计因数
      data[1] += core;
    } else {
      // 计算累计物理量
      data[0] += value * cellV * volumeF;
      // 计算累计因数
      data[1] += cellV * volumeF;
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
  static Foam::scalar rootProcOp(const base::CDTensor1& dataT) { return dataT[0]; }
};

class mixGlobalForce : public globalForce {
 public:
  cfdemTypeName("mixGlobalForce");

  cfdemDefineNewFunctionAdder(globalForce, mixGlobalForce);

  //! \brief Constructor
  mixGlobalForce(cfdemCloud& cloud);

  //! \brief Destructor
  ~mixGlobalForce();

  //! \brief 构建 expanded cell set
  void buildExpandedCellMap();

  //! \brief 每一次耦合中，在 set force 前执行
  void initBeforeSetForce();

  //! \brief 每一次耦合中，在 set force 后执行
  void endAfterSetForce();

  //! \brief 获取颗粒的 expanded cell set
  const std::unordered_set<int>& getExpandedCellSet(const int index) {
    if (expandedCellMap_.end() == expandedCellMap_.find(index)) {
      FatalError << __func__ << ": invoke buildExpandedCellMap() before get cell set" << abort(FatalError);
    }
    return expandedCellMap_[index];
  }

  //! \brief 获取颗粒处背景流体速度
  Foam::vector getBackgroundUfluid(const int index) const;

  //! \brief 获取颗粒处背景空隙率
  double getBackgroundVoidFraction(const int index) const;

  //! \brief 获取颗粒处背景流体的 ddtU
  Foam::vector getBackgroundDDtU(const int index) const;

  //! \brief 获取颗粒处背景流体的涡量
  Foam::vector getBackgroundVorticity(const int index) const;

  //! \brief 获取颗粒处背景流体的压力梯度
  Foam::vector getBackgroundGradP(const int index) const;

  //! \brief 获取颗粒处背景流体的粘性应力
  Foam::vector getBackgroundDivTau(const int index) const;

#if 1
  //! \brief 高斯核函数
  double GaussCore(const Foam::vector& particlePos, const Foam::vector& cellPos, const double radius) const {
    double dist = mag(particlePos - cellPos);
    return exp(-1.0 * dist * dist / (2 * pow(2 * radius * GaussCoreEff_, 2)));
  }
#else
  //! \brief 高斯核函数
  double GaussCore(const Foam::vector& particlePos, const Foam::vector& cellPos, const double radius) const {
    double dist = mag(particlePos - cellPos);
    return exp(-1.0 * sqr(dist) / sqr(GaussCoreEff_ * 2 * radius));
  }
#endif

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
  void setBackgroundFieldValue(const FieldType& field, std::unordered_map<int, DType>& fieldValueMap) const;

 private:
  //! \brief 背景流体速度
  //! \note map: 颗粒索引 --> 背景流体速度
  std::unordered_map<int, Foam::vector> backgroundUfluidMap_;

  //! \brief 背景流体空隙率
  //! \note map: 颗粒索引 --> 背景流体空隙率
  std::unordered_map<int, double> backgroundVoidFractionMap_;

  //! \brief 背景流体的 ddtU
  //! \note map: 颗粒索引 --> 背景流体的ddtU
  std::unordered_map<int, Foam::vector> backgroundDDtUMap_;

  //! \brief 背景流体的涡量
  //! \note map: 颗粒索引 --> 背景流体的涡量
  std::unordered_map<int, Foam::vector> backgroundVorticityMap_;

  //! \brief 背景流体的压力梯度
  //! \note map: 颗粒索引 --> 背景流体的压力梯度
  std::unordered_map<int, Foam::vector> backgroundGradPMap_;

  //! \brief 背景流体的粘性应力
  //! \note map: 颗粒索引 --> 背景流体的粘性应力
  std::unordered_map<int, Foam::vector> backgroundDivTauMap_;
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
void mixGlobalForce::setBackgroundFieldValue(const FieldType& field,
                                             std::unordered_map<int, DType>& fieldValueMap) const {
  // clear map before set
  fieldValueMap.clear();
  int number = cloud_.numberOfParticles();
  // define data buffer
  typedef base::Tensor<2, BDType> bufferTType;
  typedef base::Tensor<1, BDType> dataTType;
  bufferTType bufferTensor(base::makeShape2(number, NDim + 1), 0.0);
  // calculate data
  Foam::vector particlePos = Foam::vector::zero;  // 颗粒中心坐标
  Foam::vector cellPos = Foam::vector::zero;      // 网格中心坐标
  double radius = 0.0;                            // 颗粒半径
  double gCore = 0.0;                             // 高斯核
  double cellV = 0.0;                             // 网格体积
  for (int index = 0; index < number; ++index) {
    int findExpandedCellID = cloud_.findExpandedCellIDs()[index];
    if (findExpandedCellID >= 0) {
      particlePos = cloud_.getPosition(index);
      radius = cloud_.getRadius(index);
      auto iter = expandedCellMap_.find(index);
      if (expandedCellMap_.end() != iter) {
        // 获取颗粒扩展网格的集合
        const std::unordered_set<int>& cellIDSet = iter->second;
        for (const int& cellID : cellIDSet) {
          cellPos = cloud_.mesh().C()[cellID];
          cellV = cloud_.mesh().V()[cellID];
          // 计算高斯核
          gCore = GaussCore(particlePos, cellPos, radius);
          // 计算累计数据
          fieldRefine<DType, dataTType>::template op<fieldIsNotVF>(bufferTensor[index], field[cellID], cellV, gCore,
                                                                   voidFraction_[cellID], volumeFraction_[cellID]);
        }
      }
    }
  }
  const int masterId = 0;               // 主节点编号
  const int procId = base::procId();    // 处理器编号
  const int numProc = base::numProc();  // 处理器数量
  int tag = 100;
  // 子节点发送 bufferTensor 给主节点(非阻塞)
  if (masterId != procId) {
    MPI_Request request;
    MPI_Status status;
    base::MPI_Isend(bufferTensor, masterId, tag, &request);
    MPI_Wait(&request, &status);
  }
  // 主节点接收其他节点的 bufferTensor，并计算数据
  if (masterId == procId) {
    std::vector<MPI_Request> rVec(numProc);
    std::vector<MPI_Status> sVec(numProc);
    std::vector<bufferTType> bufferTensorVec;
    for (int inode = 0; inode < numProc; ++inode) {
      bufferTensorVec.emplace_back(base::makeShape2(number, NDim + 1), 0.0);
      if (masterId == inode) {
        base::copyTensor(bufferTensor, bufferTensorVec[inode]);
      } else {
        base::MPI_Irecv(bufferTensorVec[inode], inode, tag, rVec.data() + inode);
      }
    }
    // 主节点等待 Irecv 执行完成
    MPI_Waitall(numProc - 1, rVec.data() + 1, sVec.data() + 1);
    // 主节点重置 bufferTensor，在计算前必须重置！
    std::fill_n(bufferTensor.ptr(), bufferTensor.mSize(), 0.0);
    // 主节点计算背景流场数据
    for (int index = 0; index < number; ++index) {
      // 如果是 middle 颗粒，则计算相应的背景物理量场
      if (cloud_.checkMiddleParticle(index)) {
        for (int inode = 0; inode < numProc; ++inode) {
#pragma unroll
          for (int i = 0; i < NDim + 1; ++i) {
            bufferTensor[index][i] += bufferTensorVec[inode][index][i];
          }
        }
#pragma unroll
        for (int i = 0; i < NDim; ++i) {
          bufferTensor[index][i] /= bufferTensor[index][NDim];
        }
      }  // check middle particle
    }    // particle loop
  }
  // 主节点主节点广播 bufferTensor
  base::MPI_Bcast(bufferTensor, masterId);
  base::MPI_Barrier();
  // 计算背景流场数据
  for (int index = 0; index < number; ++index) {
    const int rootProcId = cloud_.particleRootProcIDs()[index];  // 颗粒所在的处理器编号
    if (rootProcId == procId) {
      fieldValueMap.insert(std::make_pair(index, fieldRefine<DType, dataTType>::rootProcOp(bufferTensor[index])));
    } else {
      fieldValueMap.insert(std::make_pair(index, fieldRefine<DType, dataTType>::subProcOp()));
    }
  }
  base::MPI_Barrier();
}

}  // namespace Foam

#endif  // __MIX_GLOBAL_FORCE_H__
