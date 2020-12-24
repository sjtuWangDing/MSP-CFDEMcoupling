/*!
 * Copyright (c) 2020 by Contributors
 * @file tmp.h
 * \brief Template for managing temporary objects
 *
 * @author Wang Ding
 */

#ifndef __TMP_H__
#define __TMP_H__

#include <string>
#include <type_traits>
#include <typeinfo>
#include "base/logging.h"
#include "ref_counter.h"

#ifndef TMP_DEBUG
#define TMP_DEBUG 1
#endif // TMP_DEBUG

namespace base {

//! \brief Tmp types
enum TmpType {
  //! \brief Tmp is constructed from a pointer to an object as it would be with new keyword.
  TMP = 0,
  //! \brief Tmp is constructed from a reference as would happen if the object already exists or Tmp is pushed onto
  //! stack
  CONST_REF
};

/*!
 * \brief Tmp模板的核心功能就是执行RVO(Return value optimize)，同时它还可以快速清理对象，减少峰值内存
 * Eg:
 * // Not use Tmp
 * X func() {
 *   X bigData;
 *   return bigData;
 * }
 * // Use Tmp
 * Tmp<X> func() {
 *   Tmp<X> bigData(new X);
 *   return bigData;
 * }
 * // clear Tmp
 * void func() {
 *   tmp<X> bigData(new X);
 *   bigData.clear();
 * }
 * \tparam T 被封装的类型，T类型应该具有引用计数的功能
 * \tparam RCType 引用计数器类型
 * \tparam 匿名类型形参，用于保证T类型是RCType的子类
 * \note 使用Tmp模板不一定是必须的，但是使用Tmp模板能够提高程序运行效率
 */
template <typename T, typename RCType = base::RefCounter,
          typename = typename std::enable_if<std::is_base_of<RCType, T>::value>::type>
class Tmp {
 public:
  typedef T Type;
  typedef base::RefCounter RefCounter;
  //! \brief Constructed from pointer
  inline explicit Tmp(T* ptr = 0) : type_(TMP), ptr_(ptr) {
#if TMP_DEBUG
    std::cout << "Tmp(T* ptr = 0): " << (void*)ptr_ << std::endl;
#endif
    if (ptr && !ptr->unique()) {
      CHECK(false) << __func__ << ": Attempt to construct object of type " << typeName() << " from non-unique pointer"
                   << std::endl;
    }
  }
  //! \brief Constructed from const reference
  inline Tmp(const T& ref) : type_(CONST_REF), ptr_(const_cast<T*>(&ref)) {}
  //! \brief Construct copy and increment reference count
  inline Tmp(const Tmp<T>& tmpRef) : type_(tmpRef.type_), ptr_(tmpRef.ptr_) {
#if TMP_DEBUG
    std::cout << "Tmp(const Tmp<T>& tmpRef)" << std::endl;
#endif
    if (isTmp()) {
      CHECK(isValid()) << __func__ << ": Attempt to copy an empty object of type " << typeName() << std::endl;
      operator++();
    }
  }
  //! \brief Construct copy moving content, does not increment reference count
  inline Tmp(Tmp<T>&& tmpRef) : type_(tmpRef.type_), ptr_(tmpRef.ptr_) {
#if TMP_DEBUG
    std::cout << "Tmp(Tmp<T>&& tmpRef)" << std::endl;
#endif
    if (isTmp()) {
      tmpRef.ptr_ = 0;
    }
  }
  //! \brief Construct copy transferring content of temporary if required
  inline Tmp(const Tmp<T>& tmpRef, bool allowTransfer) : type_(tmpRef.type_), ptr_(tmpRef.ptr_) {
    if (isTmp()) {
      CHECK(isValid()) << __func__ << ": Attempt to copy an empty object of type " << typeName() << std::endl;
      if (allowTransfer) {
        tmpRef.ptr_ = 0;
      } else {
        operator++();
      }
    }
  }
  //! \brief Destructor delete temporary object when the reference count is 0
  inline ~Tmp() { clear(); }
  //! \brief Assignment to pointer changing this tmp to a temporary T
  inline void operator=(T* ptr) {
    clear();
    CHECK_NE(ptr, 0) << __func__ << ": Attempt to copy an empty object of type " << typeName() << std::endl;
    CHECK(ptr->unique()) << __func__ << ": Attempt to assign of a " << typeName() << " to non-unique pointer"
                         << std::endl;
    type_ = TMP;
    ptr_ = ptr;
  }
  //! \brief Assignment transfering the temporary T to this tmp
  inline void operator=(const Tmp<T>& tmpRef) {
    clear();
    CHECK(tmpRef.isTmp()) << __func__ << ": Attempt to assign to a const reference to an object of type "
                          << typeid(T).name() << std::endl;
    CHECK(tmpRef.isValid()) << __func__ << ": Attempt to assign an empty " << typeName() << std::endl;
    type_ = TMP;
    ptr_ = tmpRef.ptr_;
    tmpRef.ptr_ = 0;
  }
  //! \brief Return true if this is really a temporary object
  inline bool isTmp() const { return TMP == type_; }
  //! \brief Return true if this temporary object empty
  inline bool isEmpty() const { return isTmp() && !ptr_; }
  //! \brief Return true if the temporary object is valid
  inline bool isValid() const { return !isTmp() || (isTmp() && ptr_); }
  //! \brief Return the type name of the TmpNRC constructed from the type name of T
  inline std::string typeName() const { return "TmpNRC<" + std::string(typeid(T).name()) + ">"; }
  //! \brief Return non-const reference or generate a fatal error if the object is const.
  inline T& ref() const {
    if (isTmp()) {
      CHECK(isValid()) << __func__ << ": Attempt to acquire reference to deallocated object with type " << typeName()
                       << std::endl;
    } else {
      CHECK(false) << __func__ << ": Attempt to acquire non-const reference to the const object with type "
                   << typeName() << std::endl;
    }
    return *ptr_;
  }
  // //! \brief Return tmp pointer for reuse.
  // inline T* ptr() {
  //   if (isTmp()) {
  //     CHECK(isValid()) << __func__ << ": Attempt to acquire reference to deallocated object with type " << typeName()
  //                      << std::endl;
  //     CHECK(ptr_->unique()) << __func__ << ": Attempt to acquire pointer to object referred to by myltiple Tmp of type "
  //                           << typeName() << std::endl;
  //     T* ptr = ptr_;
  //     ptr_ = 0;
  //     return ptr;
  //   } else {
  //     return ptr_->clone().ptr();
  //   }
  // }
  //! \brief If object pointer points to valid object: delete object and set pointer to nullptr
  inline void clear() {
    if (isTmp() && isValid()) {
      if (ptr_->unique()) {
        delete ptr_;
        ptr_ = 0;
      } else {
        ptr_->operator--();
        ptr_ = 0;
      }
    }
  }
  //! \brief Const dereference operator
  inline const T& operator()() const {
    if (isTmp()) {
      CHECK(isValid()) << __func__ << ": Attempt to acquire reference to deallocated object with type " << typeName()
                       << std::endl;
      return *ptr_;
    }
  }
  //! \brief Const cast to the underlying type reference
  inline operator const T&() const { return operator()(); }
  //! \brief Return object pointer
  inline T* operator->() {
    if (isTmp()) {
      CHECK(isValid()) << __func__ << ": Attempt to acquire reference to deallocated object with type " << typeName()
                       << std::endl;
    } else {
      CHECK(false) << __func__ << ": Attempt to cast const object to non-const pointer with type " << typeName()
                   << std::endl;
    }
    return ptr_;
  }
  //! \brief Return const object pointer
  inline const T* operator->() const {
    if (isTmp()) {
      CHECK(isValid()) << __func__ << ": Attempt to acquire reference to deallocated object with type " << typeName()
                       << std::endl;
    }
    return ptr_;
  }

 private:
  //! \brief Add ref count
  inline void operator++() {
    ptr_->operator++();
    if (ptr_->count() > 1) {
      CHECK(false) << __func__ << ": Attempt to create more than 2 tmp's referring to the same object of type "
                   << typeName() << std::endl;
    }
  }
  //! \brief Type of object
  base::TmpType type_;
  //! \brief Pointer to object
  T* ptr_;
};

}  // namespace base

#endif  // __TMP_H__
