#ifndef __GLOBAL_FORCE_H__
#define __GLOBAL_FORCE_H__

#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include "cloud/cfdem_cloud.h"
#include "mpi.h"

namespace Foam {

template <typename DType, typename TType>
struct Gauss {
  template <bool fieldIsNotVF = true>
  static void op(const TType& data, const DType& value, double cellV, double core, double vf);
  static DType subProcOp();
  static DType rootProcOp(const std::vector<TType> dataVec);
};

template <>
struct Gauss<Foam::vector, base::CDTensor1> {
  template <bool>
  static void op(const base::CDTensor1& data, const Foam::vector& value, double cellV, double core, double vf) {
    // 计算累计物理量
    data[0] += value[0] * cellV * core * vf;
    data[1] += value[1] * cellV * core * vf;
    data[2] += value[2] * cellV * core * vf;
    // 计算累计因数
    data[3] += cellV * core * vf;
  }
  static Foam::vector subProcOp() { return Foam::vector::zero; }
  static Foam::vector rootProcOp(const std::vector<base::CDTensor1>& dataVec) {
    // 由主节点计算 refined data
    Foam::vector res = Foam::vector::zero;
    double sumCore = 0.0;
    for (const auto& data : dataVec) {
      res[0] += data[0];
      res[1] += data[1];
      res[2] += data[2];
      sumCore += data[3];
    }
    res /= sumCore;
    return res;
  }
};

template <>
struct Gauss<Foam::scalar, base::CDTensor1> {
  template <bool fieldIsNotVF = true>
  static void op(const base::CDTensor1& data, const Foam::scalar& value, double cellV, double core, double vf) {
    if (fieldIsNotVF) {
      // 计算累计物理量
      data[0] += value * cellV * core * vf;
      // 计算累计因数
      data[1] += cellV * core * vf;
    } else {
      // 计算累计物理量
      data[0] += value * cellV;
      // 计算累计因数
      data[1] += cellV;
    }
  }
  static Foam::scalar subProcOp() { return 0.0; }
  static Foam::scalar rootProcOp(const std::vector<base::CDTensor1>& dataVec) {
    // 由主节点计算 refined data
    Foam::scalar res = 0.0;
    double sumCore = 0.0;
    for (const auto& data : dataVec) {
      res += data[0];
      sumCore += data[1];
    }
    res /= sumCore;
    return res;
  }
};

class globalForce {
 public:
  globalForce(cfdemCloud& cloud);

  ~globalForce();

  void updateExpandedCellMap();

  //! \brief 高斯核函数
  static double GaussCore(const Foam::vector& particlePos, const Foam::vector& cellPos, const double radius,
                          const double gaussEff) {
    double dist = mag(particlePos - cellPos);
    return exp(-1.0 * dist * dist / (2 * pow(radius * gaussEff, 2)));
  }

  template <bool fieldIsNotVF = true, int NDim = Foam::vector::dim, typename FieldType = Foam::volVectorField,
            typename DType = Foam::vector, typename BDType = Foam::scalar>
  DType getBackgroundFieldValue(int index, const FieldType& field, const volScalarField& voidFraction) const;

  inline void resetImpParticleForce() {
    impParticleForce_ == dimensionedVector("zero", impParticleForce_.dimensions(), vector::zero);
  }

  inline void resetExpParticleForcee() {
    expParticleForce_ == dimensionedVector("zero", expParticleForce_.dimensions(), vector::zero);
  }

  inline const volVectorField& impParticleForce() const { return impParticleForce_; }

  inline const volVectorField& expParticleForce() const { return expParticleForce_; }

  inline volVectorField& impParticleForce() { return impParticleForce_; }

  inline volVectorField& expParticleForce() { return expParticleForce_; }

 private:
  cfdemCloud& cloud_;

  //! \brief map: 颗粒索引 --> 扩展网格索引
  std::unordered_map<int, std::unordered_set<int>> expandedCellMap_;

  //! \brief 颗粒隐式力的总和 [N]
  volVectorField impParticleForce_;

  //! \brief 颗粒显式力的总和 [N]
  volVectorField expParticleForce_;
};

template <bool fieldIsNotVF /*= true */, int NDim /* = Foam::vector::dim*/,
          typename FieldType /* = Foam::volVectorField */, typename DType /* = Foam::vector */,
          typename BDType /* = Foam::scalar*/
          >
DType globalForce::getBackgroundFieldValue(int index, const FieldType& field,
                                           const volScalarField& voidFraction) const {
  typedef base::Tensor<1, BDType> TType;
  // 数据集合，data[0 ~ DSzie - 1] 为累计数据和，data[DSzie] 为累计因数
  TType data(base::makeShape1(NDim + 1), 0.0);
  int findExpandedCellID = cloud_.findExpandedCellIDs()[index];
  if (findExpandedCellID >= 0) {
    Foam::vector particlePos = cloud_.getPosition(index);  // 颗粒中心坐标
    Foam::vector cellPos = Foam::vector::zero;             // 网格中心坐标
    double radius = cloud_.getRadius(index);               // 颗粒半径
    double core = 0.0;                                     // 高斯核
    double cellV = 0.0;                                    // 网格体积
    if (expandedCellMap_.end() != expandedCellMap_.find(index)) {
      const std::unordered_set<int>& set = expandedCellMap_.find(index)->second;  // 颗粒扩展网格的集合
      for (int cellID : set) {
        if (cellID >= 0) {  // cell found
          cellPos = cloud_.mesh().C()[cellID];
          cellV = cloud_.mesh().V()[cellID];
          // 计算高斯核
          core = cloud_.globalF().GaussCore(particlePos, cellPos, radius, cloud_.expandedCellScale());
          // 计算累计数据
          Gauss<DType, TType>::template op<fieldIsNotVF>(data, field[cellID], cellV, core, voidFraction[cellID]);
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
    res = Gauss<DType, TType>::subProcOp();
  }
  // 主节点
  if (rootProc == procId) {
    std::vector<MPI_Request> rVec(numProc - 1);
    std::vector<MPI_Status> sVec(numProc - 1);
    std::vector<TType> dataVec;
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
    res = Gauss<DType, TType>::rootProcOp(dataVec);
  }
  // 所有处理器同步后返回结果
  base::MPI_Barrier();
  return res;
}

}  // namespace Foam

#endif  // __GLOBAL_FORCE_H__
