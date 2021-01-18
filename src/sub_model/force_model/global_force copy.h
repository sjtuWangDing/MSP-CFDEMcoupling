#ifndef __GLOBAL_FORCE_H__
#define __GLOBAL_FORCE_H__

#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include "cloud/cfdem_cloud.h"
#include "mpi.h"
#include "sub_model/data_exchange_model/data_exchange_model.h"

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

//! \brief 全局力模型，封装 force model 的共享数据与共享函数
class globalForce {
 public:
  globalForce(cfdemCloud& cloud);

  ~globalForce();

  //! \brief 每一次耦合中，在 set force 前执行
  void initBeforeSetForce();

  //! \brief 每一次耦合中，在 set force 后执行
  void endAfterSetForce();

  //! \brief 高斯核函数
  static double GaussCore(const Foam::vector& particlePos, const Foam::vector& cellPos, const double radius,
                          const double gaussEff) {
    double dist = mag(particlePos - cellPos);
    return exp(-1.0 * dist * dist / (2 * pow(radius * gaussEff, 2)));
  }

  inline void resetImpParticleForce() {
    impParticleForce_ == dimensionedVector("zero", impParticleForce_.dimensions(), vector::zero);
  }

  inline void resetExpParticleForcee() {
    expParticleForce_ == dimensionedVector("zero", expParticleForce_.dimensions(), vector::zero);
  }

  //! \brief 获取颗粒处背景流体速度
  inline Foam::vector getBackgroundUfluid(const int index) const {
    CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
    auto iter = backgroundUfluidMap_.find(index);
    if (backgroundUfluidMap_.end() != iter) {
      return iter->second;
    }
    return Foam::vector::zero;
  }

  //! \brief 获取颗粒处背景空隙率
  inline double getBackgroundVoidFraction(const int index) const {
    CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
    auto iter = backgroundVoidFractionMap_.find(index);
    if (backgroundVoidFractionMap_.end() != iter) {
      return iter->second;
    }
    return 1.0;
  }

  //! \brief 获取颗粒处背景流体的 ddtU
  inline Foam::vector getBackgroundDDtU(const int index) const {
    CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
    auto iter = backgroundDDtUMap_.find(index);
    if (backgroundDDtUMap_.end() != iter) {
      return iter->second;
    }
    return Foam::vector::zero;
  }

  //! \brief 获取上一个耦合时间步中颗粒速度
  inline Foam::vector getPrevParticleVel(const int index) const {
    // no need to check middle
    auto iter = prevParticleVelMap_.find(index);
    if (prevParticleVelMap_.end() != iter) {
      return iter->second;
    } else if (cloud_.dataExchangeM().isFirstCouplingStep()) {
      // 如果是第一个耦合时间步，则使用初始颗粒速度
      return cloud_.getInitVelocity(index);
    }
    return Foam::vector::zero;
  }

  //! \brief 获取初始时刻的相对速度 Ufluid - Up
  inline Foam::vector getInitUr(const int index) const {
    CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
    auto iter = initUrMap_.find(index);
    if (initUrMap_.end() != iter) {
      return iter->second;
    }
    return Foam::vector::zero;
  }

  //! \brief 获取颗粒的历史 ddtUr
  inline std::vector<Foam::vector>& getDDtUrHistory(const int index) {
    CHECK(cloud_.checkMiddleParticle(index)) << __func__ << ": particle " << index << " is not middle type";
    auto iter = ddtUrHistoryMap_.find(index);
    if (ddtUrHistoryMap_.end() != iter) {
      // find particle of index's history ddtUr
      return iter->second;
    }
    // not find
    ddtUrHistoryMap_.insert(std::make_pair(index, std::vector<Foam::vector>()));
    return ddtUrHistoryMap_[index];
  }

#ifdef version21
  inline const uniformDimensionedVectorField& g() const { return g_; }
#elif defined(version16ext) || defined(version15)
  inline const dimensionedVector& g() const { return g_; }
#endif

  inline const volScalarField& rhoField() const { return rho_; }

  inline const volVectorField& U() const { return U_; }

  inline const volScalarField& voidFraction() const { return voidFraction_; }

  inline const volVectorField& impParticleForce() const { return impParticleForce_; }

  inline const volVectorField& expParticleForce() const { return expParticleForce_; }

  inline const surfaceScalarField& phi() const { return phi_; }

  inline const volVectorField& ddtU() const { return ddtU_; }

  inline volVectorField& impParticleForce() { return impParticleForce_; }

  inline volVectorField& expParticleForce() { return expParticleForce_; }

 private:
  //! \brief 更新函数
  void updateKernel();

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
  cfdemCloud& cloud_;

  //! \note subPropsDict_ should be declared in front of other members
  dictionary subPropsDict_;

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

  //! \brief 上一个时间步中的颗粒的速度
  //! \note map: 颗粒索引 --> 颗粒速度
  std::unordered_map<int, Foam::vector> prevParticleVelMap_;

  //! \brief 初始时刻的相对速度 Ufluid - Up
  //! \note map: 颗粒索引 --> 初始时刻的相对速度
  std::unordered_map<int, Foam::vector> initUrMap_;

  //! \brief 记录每一个颗粒的每一个耦合时间步的 ddtUr = d(Ufluid - Up) / dt
  //! \note map: 颗粒索引 --> 历史 ddtUr
  std::unordered_map<int, std::vector<Foam::vector>> ddtUrHistoryMap_;

  //! \brief name of the gravity field
  std::string gravityFieldName_;

  //! \brief name of the rho field
  std::string densityFieldName_;

  //! \brief 速度场名称
  std::string velFieldName_;

  //! \brief 空隙率场的名称
  std::string voidFractionFieldName_;

  //! \brief 通量场的名称
  std::string phiFieldName_;

#ifdef version21
  const uniformDimensionedVectorField& g_;
#elif defined(version16ext) || defined(version15)
  const dimensionedVector& g_;
#endif

  const volScalarField& rho_;

  const volVectorField& U_;

  const volScalarField& voidFraction_;

  const surfaceScalarField& phi_;

  volVectorField ddtU_;

  //! \brief 颗粒隐式力的总和 [N]
  volVectorField impParticleForce_;

  //! \brief 颗粒显式力的总和 [N]
  volVectorField expParticleForce_;
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
DType globalForce::getBackgroundFieldValue(int index, const FieldType& field) const {
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

#endif  // __GLOBAL_FORCE_H__
