#ifndef __TENSOR_H__
#define __TENSOR_H__

#include <initializer_list>
#include <iostream>
#include "base/logging.h"
#include "base/memory/x_alloc.h"
#include "base/tensor/expression.h"

/*! \brief default data type for tensor string */
#ifndef TENSOR_DEFAULT_DTYPE
#define TENSOR_DEFAULT_DTYPE base::d_real_t
#endif

/*! \brief default device type */
#ifndef TENSOR_DEFAULT_DEVICE
#define TENSOR_DEFAULT_DEVICE base::cpu
#endif

/*! \brief default alloc type */
#ifndef TENSOR_DEFAULT_ALLOC
#define TENSOR_DEFAULT_ALLOC(DType) base::simple_alloc<DType, base::default_alloc>
#endif

/*! \brief malloc alloc type */
#ifndef TENSOR_MALLOC_ALLOC
#define TENSOR_MALLOC_ALLOC(DType) base::simple_alloc<DType, base::malloc_alloc>
#endif

/*! \brief verbose */
#ifndef TENSOR_VERBOSE
#define TENSOR_VERBOSE 0
#endif

namespace base {

/*!
 * \brief shape of a tensor
 * \tparam dimension dimension of tensor
 */
template <int dimension>
class Shape {
 public:
  /*! \brief dimension of current shape */
  static const int kDim = dimension;
  /*! \brief dimension of current shape minus one */
  static const int kSubdim = dimension - 1;
  /*! \brief storing the dimension information */
  index_t shape_[kDim];
  /*! \brief default constructor, do nothing */
  inline Shape(void) {}
  /*! \brief constuctor */
  inline Shape(const Shape<kDim>& s) {
#pragma unroll
    for (int i = 0; i < kDim; ++i) {
      this->shape_[i] = s.shape_[i];
    }
  }
  /*! \return number of valid elements */
  inline index_t size() const {
    index_t res = this->shape_[0];
    for (int i = 1; i < kDim; ++i) {
      res *= this->shape_[i];
    }
    return res;
  }
  /*!
   * \brief get corresponding index
   * \param idx dimension index
   * \return the corresponding dimension size
   */
  inline index_t& operator[](int idx) {
    CHECK(idx >= 0 && idx < kDim);
    return shape_[idx];
  }
  inline index_t operator[](int idx) const {
    CHECK(idx >= 0 && idx < kDim);
    return shape_[idx];
  }
  /*!
   * \return whether two shape equals
   * \param s the shape to compare against
   */
  inline bool operator==(const Shape<kDim>& s) const {
#pragma unroll
    for (int i = 0; i < kDim; ++i) {
      if (s.shape_[i] != this->shape_[i]) {
        return false;
      }
    }
    return true;
  }
  /*!
   * \return whether two shape not equal
   * \param s the shape to compare against
   */
  inline bool operator!=(const Shape<kDim>& s) const { return !(*this == s); }
  /*!
   * \brief get subshape that takes off largest dimension
   * \return subshape
   */
  inline Shape<kSubdim> subShape() const {
    Shape<kSubdim> s;
#pragma unroll
    for (int i = 0; i < kSubdim; ++i) {
      s.shape_[i] = this->operator[](i + 1);
    }
    return s;
  }
  /*!
   * \brief flatten the tensor, return a 1D shape
   * \return the flat 1d shape
   */
  Shape<1> flatTo1D() const {
    Shape<1> s;
    s[0] = this->size();
    return s;
  }
  /*! \brief allow string printing of the shape */
  template <int dim>
  friend std::ostream& operator<<(std::ostream& os, const Shape<dim>& shape);
};

template <int dim>
std::ostream& operator<<(std::ostream& os, const Shape<dim>& shape) {
  os << "(";
  for (int i = 0; i < dim; ++i) {
    if (i != 0) os << ", ";
    os << shape[i];
  }
  os << ")";
  return os;
}

/*!
 * \brief construct a one dimension shape, stride will equal s0
 * \param s0 size of dimension 0
 * \return the shape construction
 */
inline Shape<1> makeShape1(index_t s0) {
  Shape<1> s;
  s[0] = s0;
  return s;
}

/*!
 * \brief construct a two dimension shape, stride will equal s0
 * \param s0 size of dimension 0
 * \param s1 size of dimension 1
 * \return the shape construction
 */
inline Shape<2> makeShape2(index_t s0, index_t s1) {
  Shape<2> s;
  s[0] = s0;
  s[1] = s1;
  return s;
}

/*!
 * \brief construct shape
 * \param li initialize list
 * \return the shape construction
 */
template <index_t dimension>
Shape<dimension> makeShape(const std::initializer_list<index_t>& li) {
  CHECK_EQ(dimension, li.size());
  Shape<dimension> temp;
  int i = 0;
  for (auto it = li.begin(); it != li.end(); ++it) {
    temp[i++] = *it;
  }
  return temp;
}

/*!
 * \brief general tensor
 * \tparam dimension dimension of the tensor
 * \tparam DType the type of elements in the tensor
 * \tparam Device which device the tensor is on
 * \tparam Alloc memory allocator
 */
template <int dimension, typename DType = TENSOR_DEFAULT_DTYPE, typename Device = TENSOR_DEFAULT_DEVICE,
          typename Alloc = TENSOR_DEFAULT_ALLOC(DType)>
class Tensor : public base::Exp<Tensor<dimension, DType, Device, Alloc>, DType> {
 public:
  /*! \brief shape type of the tensor */
  typedef Shape<dimension> TShape;
  /*! \brief whether current type lies in cpu */
  static const bool kDevCPU = Device::kDevCPU;
  /*! \brief dimension */
  static const int kDim = dimension;
  /*! \brief dimension of subtype */
  static const int kSubdim = dimension - 1;
  /*! \brief pointer to the data memory allocated by self */
  DType* dptr_;
  /*! \brief pointer to the data memory allocated by other */
  DType* optr_;
  /*! \brief shape of the tensor */
  TShape shape_;
  /*!
   * \brief storing the stride information in x dimension
   *        this is used to deal with pitch allocation in gpu or sse(align x dimension to 64bit) for efficiency
   */
  index_t stride_;

  /*! \brief default constructor */
  inline Tensor() : dptr_(nullptr), optr_(nullptr), shape_(), stride_(0) {}

  /*! \brief constructor from shape */
  inline Tensor(const TShape& shape) : dptr_(nullptr), optr_(nullptr), shape_(shape), stride_(shape[kSubdim]) {
    CHECK_GE(shape.size(), 0) << "Error shape size";
    // allocate memory
    dptr_ = Alloc::allocate(shape.size());
  }

  /*! \brief constructor from shape and initial value */
  inline Tensor(const TShape& shape, DType initVal)
      : dptr_(nullptr), optr_(nullptr), shape_(shape), stride_(shape[kSubdim]) {
    CHECK_GE(shape.size(), 0) << "Error shape size";
    // allocate memory
    dptr_ = Alloc::allocate(shape.size());
    std::fill_n(dptr_, shape.size(), initVal);
  }

  /*! \brief constructor from data pointer and shape, without stride */
  inline Tensor(DType* optr, const TShape& shape)
      : dptr_(nullptr), optr_(optr), shape_(shape), stride_(shape[kSubdim]) {}

  /*! \brief constructor from data pointer and shape */
  inline Tensor(DType* optr, const TShape& shape, index_t stride)
      : dptr_(nullptr), optr_(optr), shape_(shape), stride_(stride) {}

  /*! \brief deconstructor */
  inline ~Tensor() {
    optr_ = nullptr;
    if (dptr_) {
      Alloc::deallocate(dptr_, mSize());
    }
    dptr_ = nullptr;
  }

  /*! \brief delete copy constructor */
  inline Tensor(const Tensor<dimension, DType, Device, Alloc>& that) = delete;

  /*! \brief delete copy assignment */
  inline Tensor& operator=(const Tensor<dimension, DType, Device, Alloc>& that) = delete;

  /*! \brief move constructor */
  inline Tensor(Tensor<dimension, DType, Device, Alloc>&& that) {
#if TENSOR_VERBOSE
    std::cout << "Tensor(Tensor&&)" << std::endl;
#endif
    if (this != &that) {
      dptr_ = that.dptr_;
      optr_ = that.optr_;
      shape_ = that.shape_;
      stride_ = that.stride_;
      // release old pointer
      that.dptr_ = nullptr;
      that.optr_ = nullptr;
    }
  }

  /*! \brief implement the assignment of same type */
  inline Tensor<dimension, DType, Device, Alloc>& operator=(Tensor<dimension, DType, Device, Alloc>&& that) {
#if TENSOR_VERBOSE
    std::cout << "Tensor::operator=(Tensor&&)" << std::endl;
#endif
    if (this != &that) {
      // that 为右值引用，但是本身是左值，使用 move 将其转为右值匹配移动构造
      Tensor<dimension, DType, Device, Alloc> tmp = std::move(that);
      std::swap(tmp.dptr_, dptr_);
      std::swap(tmp.optr_, optr_);
      shape_ = tmp.shape_;
      stride_ = tmp.stride_;
    }
    return *this;
  }

  /*! \brief get valid ptr */
  inline DType* ptr() const { return dptr_ ? dptr_ : optr_; }

  /*! \brief return tensor is empty */
  inline bool isEmpty() const { return nullptr == dptr_ && nullptr == optr_; }

  /*!
   * \return memory cost of the tensor, including the aligned x dimension
   * \tparam startDim the starting dimension
   */
  template <int startDim>
  inline index_t memSize() const {
    if (nullptr == dptr_ && nullptr == optr_) {
      return 0;
    }
    index_t res = this->stride_;
#pragma unroll
    for (int i = startDim; i < kSubdim; ++i) {
      res *= this->shape_[i];
    }
    return res;
  }

  /*! \return whether the tensor's memory is continuous */
  inline bool isContinus() const { return this->shape_[kSubdim] == stride_; }

  /*! \return memory cost of the tensor, including the aligned x dimension */
  inline index_t mSize() const { return this->memSize<0>(); }

  /*!
   * \brief return size of i-th dimension, start counting from highest dimension
   * \param idx the dimension count from the highest dimensin
   * \return the size
   */
  inline index_t size(int idx) const { return shape_[idx]; }

  /*!
   * \brief get a element of dimension - 1
   * \param idx index
   * \return the result tensor
   */
  inline Tensor<kSubdim, DType, Device, Alloc> operator[](int idx) const {
    CHECK(idx >= 0 && idx < static_cast<int>(shape_[0])) << " with idx = " << idx << " and shape_0 = " << shape_[0];
    // Note: not use std::move() here, because this will force invoke move constructor
    // and diable return value optimization(RVO)
    return Tensor<kSubdim, DType, Device, Alloc>(ptr() + this->memSize<1>() * idx, shape_.subShape(), stride_);
  }
};

/*! \brief respecialized class Tensor1D, thei is due to different implementation in operator[] */
template <typename DType, typename Device, typename Alloc>
class Tensor<1, DType, Device, Alloc> : public base::Exp<Tensor<1, DType, Device, Alloc>, DType> {
 public:
  typedef Shape<1> TShape;
  static const int kDim = 1;
  static const int kSubdim = 0;
  DType* dptr_;
  DType* optr_;
  Shape<1> shape_;
  index_t stride_;

  inline Tensor() : dptr_(nullptr), optr_(nullptr), shape_(), stride_(0) {}

  inline Tensor(const Shape<1>& shape) : dptr_(nullptr), optr_(nullptr), shape_(shape), stride_(shape[0]) {
    CHECK_GE(shape.size(), 0) << "Error shape size";
    // allocate memory
    dptr_ = Alloc::allocate(shape.size());
  }

  inline Tensor(const Shape<1>& shape, DType initVal)
      : dptr_(nullptr), optr_(nullptr), shape_(shape), stride_(shape[0]) {
    CHECK_GE(shape.size(), 0) << "Error shape size";
    // allocate memory
    dptr_ = Alloc::allocate(shape.size());
    std::fill_n(dptr_, shape.size(), initVal);
  }

  inline Tensor(DType* optr, const Shape<1>& shape) : dptr_(nullptr), optr_(optr), shape_(shape), stride_(shape[0]) {}

  inline Tensor(DType* optr, const Shape<1>& shape, index_t stride)
      : dptr_(nullptr), optr_(optr), shape_(shape), stride_(stride) {}

  inline ~Tensor() {
    optr_ = nullptr;
    if (dptr_) {
      Alloc::deallocate(dptr_, mSize());
    }
    dptr_ = nullptr;
  }

  inline Tensor(const Tensor<1, DType, Device, Alloc>& that) = delete;

  inline Tensor& operator=(const Tensor<1, DType, Device, Alloc>& that) = delete;

  inline Tensor(Tensor<1, DType, Device, Alloc>&& that) {
    if (this != &that) {
      dptr_ = that.dptr_;
      optr_ = that.optr_;
      shape_ = that.shape_;
      stride_ = that.stride_;
      // release old pointer
      that.dptr_ = nullptr;
      that.optr_ = nullptr;
    }
  }

  inline Tensor& operator=(Tensor<1, DType, Device, Alloc>&& that) {
    if (this != &that) {
      // that 为右值引用，但是本身是左值，使用 move 将其转为右值匹配移动构造
      Tensor<1, DType, Device, Alloc> tmp = std::move(that);
      std::swap(tmp.dptr_, dptr_);
      std::swap(tmp.optr_, optr_);
      shape_ = tmp.shape_;
      stride_ = tmp.stride_;
    }
    return *this;
  }

  /*! \brief get valid ptr */
  inline DType* ptr() const {
    CHECK(isValid());
    return dptr_ ? dptr_ : optr_;
  }

  /*! \brief return tensor is empty */
  inline bool isEmpty() const { return nullptr == dptr_ && nullptr == optr_; }

  /*! \brief return tensor is valid */
  inline bool isValid() const {
    return (nullptr != dptr_ && nullptr == optr_) || (nullptr == dptr_ && nullptr != optr_);
  }

  // member function
  template <int startDim>
  inline index_t memSize() const {
    if (nullptr == dptr_ && nullptr == optr_) {
      return 0;
    }
    CHECK_EQ(startDim, 0);
    return stride_;
  }

  inline bool isContinus() const { return stride_ == shape_[0]; }

  inline index_t mSize() const { return memSize<0>(); }

  inline index_t size(int idx) const {
    CHECK_EQ(idx, 0);
    return shape_[idx];
  }

  inline DType& operator[](int idx) const {
    CHECK(idx >= 0 && idx < static_cast<int>(mSize())) << " with idx = " << idx << " and mSize = " << mSize();
    return ptr()[idx];
  }

  friend std::ostream& operator<<(std::ostream& os, const Tensor& tensor) {
    if (!tensor.isValid()) {
      os << "[]";
    } else if (0 == tensor.mSize()) {
      os << "[]";
    } else {
      os << "[";
      for (int i = 0; i < tensor.mSize() - 1; ++i) {
        os << tensor[i] << ", ";
      }
      os << tensor[tensor.mSize() - 1] << "]";
    }
    return os;
  }
};

//! \brief define tensor type used only in CFDEM side.
using CITensor1 = Tensor<1, int, cpu, TENSOR_DEFAULT_ALLOC(int)>;
using CITensor2 = Tensor<2, int, cpu, TENSOR_DEFAULT_ALLOC(int)>;
using CDTensor1 = Tensor<1, double, cpu, TENSOR_DEFAULT_ALLOC(double)>;
using CDTensor2 = Tensor<2, double, cpu, TENSOR_DEFAULT_ALLOC(double)>;

//! \brief define tensor type used to data exchange between CFDEM and LIGGGHTS.
using CDExTensor1 = Tensor<1, double, cpu, TENSOR_MALLOC_ALLOC(double)>;
using CDExTensor2 = Tensor<2, double, cpu, TENSOR_MALLOC_ALLOC(double)>;

template <int dimension, typename DType, typename Device, typename Alloc>
void fillTensor(const Tensor<dimension, DType, Device, Alloc>& tensor, const DType& value) {
  if (tensor.isEmpty()) {
    return;
  }
  std::fill_n(tensor.ptr(), tensor.mSize(), value);
}

template <int dimension, typename DType, typename Device, typename Alloc>
void copyTensor(const Tensor<dimension, DType, Device, Alloc>& src,
                const Tensor<dimension, DType, Device, Alloc>& dest) {
  if (dest.isEmpty()) {
    return;
  }
  CHECK_EQ(dest.mSize(), src.mSize()) << ": destination tensor's size should equal to src tensor'size";
  for (index_t i = 0; i < dest.mSize(); ++i) {
    dest.ptr()[i] = src.ptr()[i];
  }
}

}  // namespace base

#endif  // __TENSOR_H__
