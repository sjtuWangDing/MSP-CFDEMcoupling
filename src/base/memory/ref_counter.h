/*!
 * Copyright (c) 2020 by Contributors
 * @file ref_counter.h
 * \brief Reference counter
 *
 * @author Wang Ding
 */

#ifndef __REF_COUNTER_H__
#define __REF_COUNTER_H__

#include <stdint.h>

namespace base {

typedef int32_t Atomic32;

/*!
 * \brief 引用计数
 * \note 如果希望使用Tmp模板封装某一个类，那么这个类就需要继承RefCount
 */
class RefCounter {
 private:
  Atomic32 count_;
  RefCounter(const RefCounter&);
  RefCounter& operator=(const RefCounter&);

 protected:
  RefCounter(int count = 0) : count_(count) {}

 public:
  int count() const { return count_; }
  bool unique() const { return 0 == count_; }
  RefCounter& operator++() {
    __sync_fetch_and_add(&count_, 1);
    return *this;
  }
  RefCounter& operator--() {
    __sync_fetch_and_sub(&count_, 1);
    return *this;
  }
  RefCounter operator++(int) {
    RefCounter res(count_);
    this->operator++();
    return res;
  }
  RefCounter operator--(int) {
    RefCounter res(count_);
    this->operator--();
    return res;
  }
};

}  // namespace base

#endif  // __REF_COUNTER_H__
