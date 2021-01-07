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
  dataExchangeModel
\*---------------------------------------------------------------------------*/

#ifndef __DATA_EXCHANGE_MODEL_H__
#define __DATA_EXCHANGE_MODEL_H__

#include "base/run_time_selection_tables.h"
#include "cloud/cfdem_cloud.h"

namespace Foam {

class dataExchangeModel {
 public:
  //! \brief Runtime type information
  cfdemBaseTypeName("dataExchangeModel", "");

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(autoPtr, dataExchangeModel, (cfdemCloud & cloud), (cloud));

  //! \brief Selector
  static autoPtr<dataExchangeModel> New(cfdemCloud& cloud, const dictionary& dict);

  //! \brief Constructor
  dataExchangeModel(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~dataExchangeModel();

  //! \brief 耦合时间步计数器
  class CouplingStepCounter {
   public:
    CouplingStepCounter(const dataExchangeModel& model) : model_(model) {}
    ~CouplingStepCounter() { const_cast<dataExchangeModel&>(model_).couplingStep_ += 1; }
    CouplingStepCounter(const CouplingStepCounter&) = delete;
    CouplingStepCounter& operator=(const CouplingStepCounter&) = delete;
    CouplingStepCounter(CouplingStepCounter&&) = delete;
    CouplingStepCounter& operator=(CouplingStepCounter&&) = delete;

   private:
    const dataExchangeModel& model_;
  };

  /*!
   * \brief 检查时间步长是否满足耦合要求（在初始化模型的时候调用）
   * \note (1) 耦合时间步长应该大于 or 等于 CFD 时间步长
   *       (2) 耦合时间步长应该为 CFD 时间步长的整数倍
   *       (3) 如果耦合时间步长 != CFD 时间步长，则要求 allowUseSubCFDTimeStep() 为 true
   */
  void checkTimeStepSize() const;

  /*!
   * \brief 因为耦合时间步长 = 流体时间步长的整数倍，所以 timeStepFraction()
   *   用于计算每个流体时间步在耦合时间步中的所占比例，如果 couplingTime() == 3 * CFDts，
   *   那么每一个耦合时间步由 3 个流体时间步构成，那么这三个流体时间步的
   *   timeStepFraction() 分别返回 0, 0.333333, 0.666666
   */
  double timeStepFraction() const;

  /*!
   * \brief 因为耦合时间步长 = 流体时间步长的整数倍，而耦合发生在耦合时间步中的第一个流体时间步中，
   *        所以判断当前流体时间步是否同时也是耦合时间步
   */
  bool checkValidCouplingStep() const;

  template <int nDim, typename DType, typename Device, typename Alloc>
  void free(base::Tensor<nDim, DType, Device, Alloc>& tensor, DType**& pData) const {
    CHECK(2 == nDim || 1 == nDim) << "tensor free: error tensor dimension";
    // 这里不能使用二级配置器 base::default_alloc
    // 使用一级配置器 base::malloc_alloc，因为通过 LIGGGHTS 分配的内存需要一级配置器释放
    CHECK((std::is_same<Alloc, TENSOR_MALLOC_ALLOC(DType)>::value))
        << "dataExchangeModel::free: please use malloc_alloc because LIGGGHTS using malloc() to allocate memory";
    if (nullptr == pData && tensor.isEmpty()) {
      // no need to free
      return;
    } else if (nullptr != pData && nullptr != pData[0]) {
      CHECK_EQ(nullptr, tensor.dptr_)
          << "dataExchangeModel::free: tensor's dptr_ should be nullptr because of using liggghts's memory allocate";
      CHECK_EQ(pData[0], tensor.optr_) << "dataExchangeModel::free: tensor's optr_ not match with pData[0]";
      // 释放 pData[0] 指向的内存
      TENSOR_MALLOC_ALLOC(DType)::deallocate(pData[0], tensor.mSize());
      // 重置 pData[0]
      pData[0] = nullptr;
      // 释放 pData 指向的内存
      TENSOR_MALLOC_ALLOC(DType*)::deallocate(pData, tensor.size(0));
      // Note: reset pointer
      // Note: must reset pointer of pData to nullptr, otherwise liggghts will allocate memory failed!
      pData = nullptr;
      tensor.optr_ = nullptr;
    } else {
      CHECK(false) << "dataExchangeModel::free: bad free";
    }
  }

  /*!
   * \brief 重新分配 tensor 内存（当前 tensor 用于与 LIGGGHTS 交换数据，所以通过 LIGGGHTS 提供的接口分配内存）
   * \param tensor 需要分配的内存对象
   * \param shape 指定 tensor 的 shape
   * \param pData 与该 tensor 强绑定的 pointer
   * \param initVal 初始值
   */
  template <int nDim, typename DType, typename Device, typename Alloc>
  void realloc(base::Tensor<nDim, DType, Device, Alloc>& tensor, const base::Shape<nDim>& shape, DType**& pData,
               DType initVal) const {
    CHECK(2 == nDim || 1 == nDim) << "tensor realloc: error tensor dimension";
    this->free(tensor, pData);
    index_t nrow = shape[0];
    index_t ncol = 1 == nDim ? 1 : shape[1];
    // 调用 LIGGGHTS 分配内存到 pData
    liggghtsAllocate(pData, nrow, ncol, initVal);
    // 将 pData[0] 转移至 tensor
    // 这里会匹配 tensor 中的 move assignment
    tensor = std::move(base::Tensor<nDim, DType, Device, Alloc>(pData[0], shape));
  }

  //! \return 当前耦合时间步中颗粒的数量
  virtual int couple() = 0;

  /*!
   * \brief 从 LAMMPS 中获取数据（数据类型为 double）
   *        Eg: getData("x", "vector-atom", array);
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  virtual void getData(const std::string& dataName, const std::string& dataType, double** array) const = 0;

  /*!
   * \brief 从 LAMMPS 中获取数据（数据类型为 int）
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  virtual void getData(const std::string& dataName, const std::string& dataType, int** array) const = 0;

  /*!
   * \brief 传递数据到 LAMMPS 中（数据类型为 double）
   *        Eg: giveData("dragforce", "vector-atom", array);
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  virtual void giveData(const std::string& dataName, const std::string& dataType, double** field) const = 0;

  /*!
   * \brief 传递数据到 LAMMPS 中（数据类型为 int）
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  virtual void giveData(const std::string& dataName, const std::string& dataType, int** field) const = 0;

  //! \brief Allocate for 2-D double array using liggghts interface
  virtual void liggghtsAllocate(double**& array, int length, int width, double initVal = 0.0) const = 0;

  //! \brief Allocate 2-D int array using liggghts interface
  virtual void liggghtsAllocate(int**& array, int length, int width, int initVal = 0) const = 0;

  // //! \brief Allocate and destroy for 2-D double array
  // virtual void destroy(double**& array) {};
  // virtual void destroy(double**& array, int) {};
  // virtual void allocateArray(double**& array, double initVal, int width, int length) {};
  // virtual void allocateArray(double**& array, double initVal, int width, const char* length = "nparticles") {};

  // //! \brief Allocate and destroy for 2-D int array
  // virtual void destroy(int**& array) {};
  // virtual void destroy(int**& array, int) {};
  // virtual void allocateArray(int**& array, int initVal, int width, int length) {};
  // virtual void allocateArray(int**& array, int initVal, int width, const char* length  = "nparticles") {};

  // //! \brief Allocate and destroy for 1-D double array
  // virtual void destroy(double*& array) {};
  // virtual void allocateArray(double*& array, double initVal, int length) {};

  // //! \brief Allocate and destroy for 1-D int array
  // virtual void destroy(int*& array) {};
  // virtual void allocateArray(int*& array, int initVal, int length) {};

  // virtual void getData(const std::string& dataName,
  //                      const std::string& dataType,
  //                      double** const& field,
  //                      label step) = 0;
  // virtual void getData(const std::string& dataName,
  //                      const std::string& dataType,
  //                      int** const& field,
  //                      label step) = 0;
  // virtual void giveData(const std::string& dataName,
  //                       const std::string& dataType,
  //                       double** const& field,
  //                       const char* fieldType = "double") = 0;

  // //! \brief 释放离散内存
  // virtual void destroyDiscreteMemory(double** const& array, int length) {
  //   FatalError << "dataExchangeModel::destroyDiscreteMemory(): "
  //     << "using base class function, please use derived class function\n"
  //     << abort(FatalError);
  // }
  // //! \brief 释放离散内存
  // virtual void destroyDiscreteMemory(int** const& array, int length) {
  //   FatalError << "dataExchangeModel::destroyDiscreteMemory(): "
  //     << "using base class function, please use derived class function\n"
  //     << abort(FatalError);
  // }

  inline int couplingStep() const { return couplingStep_; }

  inline double DEMts() const { return DEMts_; }

  //! \brief 计算耦合间隔的秒数
  inline double couplingTime() const { return cloud_.cProps().couplingInterval() * DEMts_; }

  //! \brief 当前耦合时间步的起始时间（全局时间）
  inline double TSstart() const {
    return cloud_.mesh().time().startTime().value() + (couplingStep_ - 1) * couplingTime();
  }

  //! \brief 当前耦合时间步的结束时间（全局时间）
  inline double TSend() const { return cloud_.mesh().time().startTime().value() + couplingStep_ * couplingTime(); }

  //! \brief 计算指定时间的 DEM 步
  inline int DEMstepsTillT(double t) const {
    return (t - (cloud_.mesh().time().value() - couplingTime()) + SMALL) / DEMts_;
  }

 protected:
  cfdemCloud& cloud_;

  //! \brief 初始化 dataExchangeModel 的时候，即耦合开始的时候，记录当前的流体时间步，通常从 0 开始
  const int timeIndexOffset_;

  //! \brief 耦合时间步的计数，每执行一次 coupling 记录一次
  int couplingStep_;

  //! \brief DEM 的时间步长
  scalar DEMts_;
};

}  // namespace Foam

#endif  // __DATA_EXCHANGE_MODEL_H__