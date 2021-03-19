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
#include <cstdlib>
#include <string.h>

namespace base {

//! \brief 字符串修剪，截去字符串的首尾空白字符(空格、制表、回车、换行等)
inline std::string& trim(std::string& str) {
  // 在字符串中查找第一个与" \t\r\n"都不匹配的字符，返回它的位置，如果没有找到，返回 string::npos
  std::string::size_type first = str.find_first_not_of(" \t\r\n");
  // 在字符串中查找最后一个与" \t\r\n"都不匹配的字符，返回它的位置，如果没有找到，返回 string::npos
  std::string::size_type last = str.find_last_not_of(" \t\r\n");

  if (std::string::npos == first || std::string::npos == last) {
    // 如果没有找到任何非空白字符串，说明 str 要么为空，要么全为空白字符，则直接清空
    str.clear();
  } else {
    // 截取从第一个到最后一个非空白字符之间的子串
    str = str.substr(first, last - first + 1);
  }
  // 返回修改后的字符串
  return str;
}

/*!
 * \brief 字符串拆分，以delim字符串中的字符作为分隔符，对str字符串进行拆分，
 *   并对每个被拆分出的子串做修剪，拆分次数不超过limit，除非该参数的值为0
 * \return 被拆分出的子字符串
 */
std::vector<std::string> split(const std::string& str, const std::string& delim, int limit = 0) {
  // 存放拆分结果的子字符串
  std::vector<std::string> strv;

  // 多分配一个字符用于存放'\0'
  char temp[str.size() + 1];
  // strtok()要求待拆分字符串可写，所以将字符串拷贝到可写内存 temp 中
  strcpy(temp, str.c_str());

  // delim : " ,.:;()$"
  // limit : 3
  // temp : The quick,brown.fox:jumps;over(the)lazy$dog/0
  // strv : The
  //            quick
  //                  brown
  //                        fox:jumps;over(the)lazy$dog
  // --limit : 2   1     0

  // 依次提取子字符串
  // 第一次调用 strtok(temp, delim.c_str())，strtok函数开始遍历 temp，一旦出现 delim.c_str()
  // 中出现的字符，strtok函数停止，将 temp 中的那个字符变为'\0'，返回开始搜索的指针 temp
  // 第二次调用 strtok(NULL, delim.c_str())，则strtok开始从变为'\0'字符的下一个字符开始寻找第一个落在 delim.c_str()
  // 中的字符
  // 当最后一个子串被提取出来后，再次调用strtok是发现是从字符串的结尾'\0'开始寻找的，这说明字符串已经找完了，则返回NULL
  for (char* token = strtok(temp, delim.c_str()); token != NULL; token = strtok(NULL, delim.c_str())) {
    // 拆分后的子串
    std::string part(token);
    // 修剪后存入子串向量中
    strv.push_back(trim(part));
    // 如果拆分次数已到(若limit取缺省值0，则if条件永远为false)
    // token += strlen(token) 此时 token 位于下一个子串的开始位置
    // token - temp 返回已经拆分的字符串长度，如果小于总长度，则说明本次拆出的子串不是最后一个子串
    if (!(--limit) && (token += strlen(token)) - temp < (int)str.size()) {
      // 此时 token 位于本次拆出的子串后的'\0'位置
      // 将 ++token 以及其后的所有字符赋值到part中
      part = ++token;
      strv.push_back(trim(part));
      // 提前结束循环
      break;
    }
  }
  // 返回拆分结果的子串向量
  return strv;
}

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

}  // namespace base

#endif  // __STRING_UTILS_H__
