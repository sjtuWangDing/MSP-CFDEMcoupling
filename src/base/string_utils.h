/*!
 * Copyright (c) 2020 by Contributors
 * @file string_utils.h
 * \brief Basic functions of string utils.
 *
 * @author Wang Ding
 */

#ifndef __STRING_UTILS_H__
#define __STRING_UTILS_H__

#include <sstream>
#include <string>

inline void makeStringKernel(std::stringstream& ss) {}

template <typename DType>
inline void makeStringKernel(std::stringstream& ss, const DType& data) {
  ss << data;
}

template <typename DType, typename... ATypes>
inline void makeStringKernel(std::stringstream& ss, const DType& data, const ATypes&... args) {
  makeStringKernel(ss, data);
  makeStringKernel(ss, args...);
}

template <typename... ATypes>
std::string makeString(const ATypes&... args) {
  std::stringstream ss;
  makeStringKernel(ss, args...);
  return std::string(ss.str());
}

//! \brief Specializations for std::string type.
template <>
inline std::string makeString(const std::string& str) {
  return str;
}

inline std::string makeString(const char* cstr) {
  return std::string(cstr);
}

#endif  // __STRING_UTILS_H__
