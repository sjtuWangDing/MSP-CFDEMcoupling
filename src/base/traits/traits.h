/*!
 * Copyright (c) 2020 by Contributors
 * @file traits.h
 * \brief Traits class for primitives.
 * \note All primitives need a specialised version of this class. The
         specialised version will normally also require a conversion method.
 *
 * @author Wang Ding
 */

#ifndef __TRAITS_H__
#define __TRAITS_H__

#include <iostream>

namespace base {

template <typename PrimitiveType>
class Traits : public PrimitiveType {
 public:
  //! \brief Construct from primitive type
  Traits(const PrimitiveType& p) : PrimitiveType(p) {}
  //! \brief Construct from Istream
  Traits(std::istream& is) : PrimitiveType(is) {}
};

}  // namespace base

#endif  // __TRAITS_H__
