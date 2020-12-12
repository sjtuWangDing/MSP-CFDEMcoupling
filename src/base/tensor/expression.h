#ifndef __EXPRESSION_H__
#define __EXPRESSION_H__

namespace base {

/*!
 * \brief base class for expression
 * \tparam SubType inheritated class must put their type into this parameter
 * \tparam DType the data type of each element in the expression
 */
template <typename SubType, typename DType>
struct Exp {
 public:
  /*! \return  subtype instance of current class */
  inline const SubType& self(void) const { return *static_cast<const SubType*>(this); }
  /*! \return reference of subtype instance of current class */
  inline SubType* ptrself(void) { return static_cast<SubType*>(this); }
};

}  // namespace base

#endif  // __EXPRESSION_H__
