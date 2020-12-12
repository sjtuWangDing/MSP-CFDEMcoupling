#ifndef __TYPEINFO_H__
#define __TYPEINFO_H__

#include <stdexcept>
#include <typeinfo>
#include "base/logging.h"

namespace base {

//! \brief Type cast template function
//    wraps dynamic_cast to handle bad_cast exception.
template <typename ToType, typename FromType>
inline ToType dynamicCast(FromType from) {
  try {
    return dynamic_cast<ToType>(from);
  } catch (const std::bad_cast& ex) {
    CHECK(false) << "Attempt to cast type " << typeid(FromType).name() << " to type " << typeid(ToType).name();
  }
  return dynamic_cast<ToType>(from);
}

//! \brief Reference type cast template function
//    wraps dynamic_cast to handle bad_cast exception,
//    but handles type names via the virtual typeName() method.
template <typename ToType, typename FromType>
inline ToType& refCast(FromType& from) {
  try {
    return dynamic_cast<ToType&>(from);
  } catch (const std::bad_cast& ex) {
    CHECK(false) << "Attempt to cast type " << from.typeName() << " to type " << ToType::cTypeName();
  }
  return dynamic_cast<ToType&>(from);
}

//! \brief Pointer type cast template function
//    wraps dynamic_cast to handle bad_cast exception,
//    but handles type names via the virtual type() method.
template <typename ToType, typename FromType>
inline ToType* refCast(FromType* from) {
  try {
    return dynamic_cast<ToType*>(from);
  } catch (const std::bad_cast& ex) {
    CHECK(false) << "Attempt to cast type " << from.typeName() << " to type " << ToType::cTypeName();
  }
  return dynamic_cast<ToType*>(from);
}

}  // namespace base

#endif  // __TYPEINFO_H__
