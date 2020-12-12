#ifndef __X_ALLOC_H__
#define __X_ALLOC_H__

#include <cstdlib>

#ifndef __THROW_BAD_ALLOC
#include <new>
#define __THROW_BAD_ALLOC throw std::bad_alloc()
#endif  // __THROW_BAD_ALLOC

// window thread
#ifdef __WIN32_THREADS
#include <windows.h>
#define __NODE_ALLOCATOR_LOCK \
  if (__threads) EnterCriticalSection(&_S_node_allocator_lock)
#define __NODE_ALLOCATOR_UNLOCK \
  if (__threads) LeaveCriticalSection(&_S_node_allocator_lock)
#define __NODE_ALLOCATOR_THREADS true
#endif  // __WIN32THREADS

// posix thread
#ifdef __P_THREADS
#include <pthread.h>
#define __NODE_ALLOCATOR_LOCK \
  if (__threads) pthread_mutex_lock(&_S_node_allocator_lock)
#define __NODE_ALLOCATOR_UNLOCK \
  if (__threads) pthread_mutex_unlock(&_S_node_allocator_lock)
#define __NODE_ALLOCATOR_THREADS true
#endif  // __PTHREADS

// c++11 stl thread
#ifdef __STL_THREADS
#include <mutex>
#define __NODE_ALLOCATOR_LOCK \
  if (__threads) _S_node_allocator_lock.lock();
#define __NODE_ALLOCATOR_UNLOCK \
  if (__threads) _S_node_allocator_lock.unlock();
#define __NODE_ALLOCATOR_THREADS true
#endif  // __STL_THREADS

#if !defined(__NO_THREADS) && !defined(__WIN32_THREADS) && !defined(__P_THREADS) && !defined(__STL_THREADS)
#define __NO_THREADS
#endif

#ifdef __NO_THREADS
#define __NODE_ALLOCATOR_LOCK
#define __NODE_ALLOCATOR_UNLOCK
#define __NODE_ALLOCATOR_THREADS false
#endif

namespace base {

/*!
 * \brief 内存配置器接口（提供给外部调用）
 * \tparam _T 分配内存的类型
 * \tparam _Alloc 指定内存配置器，可以是一级或者二级配置器
 */
template <typename _T, typename _Alloc>
class simple_alloc {
 public:
  static _T* allocate(size_t _N) { return (_T*)(0 == _N ? 0 : _Alloc::allocate(_N * sizeof(_T))); }
  static _T* allocate() { return (_T*)_Alloc::allocate(sizeof(_T)); }
  static void deallocate(_T* _p, size_t _N) {
    if (0 != _N) _Alloc::deallocate(_p, _N * sizeof(_T));
  }
  static void deallocate(_T* _p) { _Alloc::deallocate(_p, sizeof(_T)); }
};

//! \brief 定义内存分配失败后 handler 的类型
typedef void (*__oom_handler)();

/*!
 * \brief SGI STL 一级配置器
 * \tparam __inst 目前暂时没有作用
 */
template <int __inst>
class __malloc_alloc_template {
 private:
  static __oom_handler _S_oom_handler;
  // 以下函数将用来处理内存不足的情况
  static void* _S_oom_malloc(size_t);
  static void* _S_oom_realloc(void*, size_t);

 public:
  //! \brief 一级配置器直接调用 malloc 分配内存
  static void* allocate(size_t _N) {
    void* _result = malloc(_N);
    // 当内存分配失败的时候，调用 _S_oom_malloc
    if (0 == _result) {
      _result = _S_oom_malloc(_N);
    }
    return _result;
  }
  //! \brief 一级配置器直接调用 realloc 分配内存（这里多了一个哑元参数是为了与二级配置器中 reallocate 函数形式统一）
  static void* reallocate(void* _p, size_t /* _old_size */, size_t _new_size) {
    void* _result = realloc(_p, _new_size);
    // 当内存分配失败的时候，调用 _S_oom_realloc
    if (0 == _result) {
      _result = _S_oom_realloc(_p, _new_size);
    }
    return _result;
  }
  //! \brief 一级配置器直接调用 free 释放内存
  static void deallocate(void* _p, size_t /* _N */) { free(_p); }
  /*!
   * \brief This static function is defined to simulate std::set_new_handler() in C++
   * \note Why not use std::set_new_handler()? 因为一级配置器并没有使用 ::operator new() 来分配内存
   */
  static void (*set_malloc_handler(void (*_f)()))() {
    void (*_old_handler)() = _S_oom_handler;
    _S_oom_handler = _f;
    return _old_handler;
  }
};

template <int __inst>
__oom_handler __malloc_alloc_template<__inst>::_S_oom_handler = 0;

template <int __inst>
void* __malloc_alloc_template<__inst>::_S_oom_malloc(size_t _N) {
  __oom_handler _my_handler = 0;
  void* _result = 0;
  // 不断尝试释放内存，分配内存
  for (;;) {
    _my_handler = _S_oom_handler;
    if (0 == _my_handler) {
      __THROW_BAD_ALLOC;
    }
    // invoke handler, try to dealloc memory
    _my_handler();
    // try to malloc memory
    _result = malloc(_N);
    if (_result) {
      return _result;
    }
  }
}

/*!
 * \brief realloc memory on pointer
 * \param _p memory address for realloc memory
 * \param _N realloc memory's size
 */
template <int __inst>
void* __malloc_alloc_template<__inst>::_S_oom_realloc(void* _p, size_t _N) {
  __oom_handler _my_handler = 0;
  void* _result = 0;
  for (;;) {
    _my_handler = _S_oom_handler;
    if (0 == _my_handler) {
      __THROW_BAD_ALLOC;
    }
    // invoke handler, try to dealloc memory
    _my_handler();
    // try to realloc memory
    _result = realloc(_p, _N);
    if (_result) {
      return _result;
    }
  }
}

// 直接将参数 __inst 指定为0
typedef __malloc_alloc_template<0> malloc_alloc;

/*!
 * \brief SGI STL 二级配置器，GCC 默认使用的是二级配置器，使用二级配置器可以避免小的内存申请造成内存碎片
 * \tparam __threads 多线程环境使用
 * \tparam __inst 目前暂时没有作用
 */
template <bool __threads, int __inst>
class __default_alloc_template {
 public:
  /*!
   * \brief 二级配置器申请内存函数
   * \note 该函数用于外部调用
   * \param _N 申请内存大小
   * \return 内存地址
   */
  static void* allocate(size_t _N) {
    void* _result = 0;
    // 如果内存块的需求大于 128 bytes，则直接调用一级配置器分配内存
    if (_N > _S_max_bytes) {
      _result = base::malloc_alloc::allocate(_N);
    } else {
      // 在 16 个 free list 中寻找对应的空闲链表
      __Obj* volatile* _my_free_list = _S_free_list + _S_free_list_index(_N);
#ifndef __NO_THREADS
      // 如果支持多线程则定义 lock
      __Lock _lock_instance;
#endif  // __NO_THREADS
      __Obj* _new_node = *_my_free_list;
      if (0 == _new_node) {
        // free list 没有可用数据块，就将区块大小先调整至 8 倍数边界，然后调用 _S_refill() 重新填充
        _result = _S_refill(_S_round_up(_N));
      } else {
        // 如果 free list 中有空闲数据块，则取出一个，并把空闲链表的指针指向下一个数据块
        *_my_free_list = _new_node->_next;
        _result = _new_node;
      }
    }
    return _result;
  }
  /*!
   * \brief 二级配置器释放内存函数
   * \note 该函数用于外部调用
   * \param _p 内存地址
   * \param _N 内存大小
   */
  static void deallocate(void* _p, size_t _N) {
    if (_N > _S_max_bytes) {
      // 大于 128 bytes，则调用第一级配置器释放
      base::malloc_alloc::deallocate(_p, _N);
    } else {
      // 将空间回收到相应空闲链表中
      __Obj* volatile* _my_free_list = _S_free_list + _S_free_list_index(_N);
#ifndef __NO_THREADS
      __Lock _lock_instance;
#endif  // __NO_THREADS
      // 将 _p 指向的空闲内存插入到相应空闲链表的表头
      // 因为 _p 指向的空闲内存的大小一定是 8 的倍数，所以利用空闲内存的前 8 个字节存放 *_my_free_list
      ((__Obj*)_p)->_next = *_my_free_list;
      *_my_free_list = (__Obj*)_p;
    }
  }

 private:
  /*!
   * \brief 自由链表中的节点结构
   *        Eg:
   *        Obj* free_list = 0;         // 链表头
   *        void* memory1 = malloc(32); // 第一个节点中的内存
   *        void* memory2 = malloc(32); // 第二个节点中的内存
   *        void* memory3 = malloc(32); // 第三个节点中的内存
   *
   *        // 将 memory1 memory2 memory3 构成链表的本质就是：
   *        // 使用 memory1 指向空间的前 8 个字节存储下一个节点的地址
   *        free_list = (Obj*)memory1;
   *        ((Obj*)memory1)->next = (Obj*)memory2;
   *        ((Obj*)memory2)->next = (Obj*)memory3;
   *        ((Obj*)memory3)->next = (Obj*)0;
   *
   *        // 获取 free_list 链表中的每个节点的内存
   *        Obj* obj1 = free_list;
   *        Obj* obj2 = obj1->next;
   *        Obj* obj3 = obj2->next;
   *        // data 字符数组中第一个字符就位于 obj1 中第一个字节，所以 data 就是 obj1 中内存的地址
   *        char* p1 = obj1->data;
   *        char* p2 = obj2->data;
   *        char* p3 = obj3->data;
   *        // 或者可以直接转换指针
   *        void* p1 = (void*)obj1;
   *        void* p2 = (void*)obj2;
   *        void* p3 = (void*)obj3;
   */
  union __Obj {
    union __Obj* _next;
    char _data[1];
  };
  //! \brief 自由链表以 8 字节为对齐方式，一直扩充到 128
  enum { _ALIGN_BYTES = 8 };
  //! \brief 内存池最大的 chunk
  enum { _MAX_BYTES = 128 };
  //! \brief 自由链表的个数
  enum { _FREE_LIST_NUM = _MAX_BYTES / _ALIGN_BYTES };
  //! \brief size_t 类型的 _ALIGN_BYTES
  static size_t _S_align_bytes;
  //! \brief size_t 类型的 _MAX_BYTES
  static size_t _S_max_bytes;
  //! \brief size_t 类型的 _FREE_LIST_NUM
  static size_t _S_free_list_num;
  /*!
   * \brief 维护 16 个空闲链表 free-list，初始化为 0，即每个链表中都没有空闲数据块
   *        各自管理大小分别为 8， 16， 24， 32，...128 bytes的内存块
   */
  static __Obj* volatile _S_free_list[];
  //! \brief 内存池中空闲内存的起始位置
  static char* _S_free_start;
  //! \brief 内存池中空闲内存的结束位置
  static char* _S_free_end;
  //! \brief 内存池的空闲内存大小
  static size_t _S_heap_size;

//! \brief 互斥量，保证多线程情况下访问 _S_free_list 的安全
#ifdef __P_THREADS
  static pthread_mutex_t _S_node_allocator_lock;
#endif  // __P_THREADS

#ifdef __STL_THREADS
  static std::mutex _S_node_allocator_lock;
#endif  // __STL_THREADS

#ifdef __WIN32_THREADS
  static CRITICAL_SECTION _S_node_allocator_lock;
  static bool _S_node_allocator_lock_initialized;

 public:
  //! \brief 只有在 __WIN32_THREADS 下需要添加构造函数
  __default_alloc_template() {
    if (false == _S_node_allocator_lock_initialized) {
      InitializedCriticalSection(&_S_node_allocator_lock);
      _S_node_allocator_lock_initialized = true;
    }
  }

 private:
#endif  // __WIN32_THREADS

  //! \brief lock used for __default_alloc_template
  class __Lock {
   public:
    __Lock() { __NODE_ALLOCATOR_LOCK; }
    ~__Lock() { __NODE_ALLOCATOR_UNLOCK; }
  };
  friend class __Lock;

  //! \brief 将任何小的内存需求上调至 8 的倍数（8, 16, 24, ... 128）
  static inline size_t _S_round_up(size_t _bytes) { return (_bytes + _S_align_bytes - 1) & (~(_S_align_bytes - 1)); }
  //! \brief 根据申请内存块的大小找到对应空闲链表的下标
  static inline size_t _S_free_list_index(size_t _bytes) { return (_bytes + _S_align_bytes - 1) / _S_align_bytes - 1; }
  /*!
   * \brief 申请 _nobjs 个 _size 大小的节点，同时最终节点数量可能被减少
   * \param _size 申请的节点大小
   * \param _nobjs [in, out] 最终申请到的节点个数
   */
  static char* _S_chunk_alloc(size_t _size, int& _nobjs);
  /*!
   * \brief 如果管理大小为 _N 的 free list 中没有可用的节点，则该函数会为该 free list 填充多个新的节点
   * \param _N 指定 free list 管理的节点大小
   */
  static void* _S_refill(size_t _N);
};

template <bool __threads, int __inst>
size_t __default_alloc_template<__threads, __inst>::_S_align_bytes = (size_t)
    __default_alloc_template<__threads, __inst>::_ALIGN_BYTES;

template <bool __threads, int __inst>
size_t __default_alloc_template<__threads, __inst>::_S_max_bytes = (size_t)
    __default_alloc_template<__threads, __inst>::_MAX_BYTES;

template <bool __threads, int __inst>
size_t __default_alloc_template<__threads, __inst>::_S_free_list_num = (size_t)
    __default_alloc_template<__threads, __inst>::_FREE_LIST_NUM;

template <bool __threads, int __inst>
typename __default_alloc_template<__threads, __inst>::__Obj* volatile __default_alloc_template<
    __threads, __inst>::_S_free_list[_FREE_LIST_NUM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

template <bool __threads, int __inst>
char* __default_alloc_template<__threads, __inst>::_S_free_start = 0;

template <bool __threads, int __inst>
char* __default_alloc_template<__threads, __inst>::_S_free_end = 0;

template <bool __threads, int __inst>
size_t __default_alloc_template<__threads, __inst>::_S_heap_size = 0;

#ifdef __P_THREADS
template <bool __threads, int __inst>
pthread_mutex_t __default_alloc_template<__threads, __inst>::_S_node_allocator_lock = PTHREAD_MUTEX_INITIALIZER;
#endif  // __P_THREADS

#ifdef __STL_THREADS
template <bool __threads, int __inst>
std::mutex __default_alloc_template<__threads, __inst>::_S_node_allocator_lock;
#endif  // __STL_THREADS

#ifdef __WIN32_THREADS
template <bool __threads, int __inst>
CRITICAL_SECTION __default_alloc_template<__threads, __inst>::_S_node_allocator_lock;

template <bool __threads, int __inst>
bool __default_alloc_template<__threads, __inst>::_S_node_allocator_lock_initialized = false;
#endif  // __WIN32_THREADS

/*!
 * \brief 申请 _nobjs 个 _size 大小的节点，同时最终节点数量可能被减少
 * \param _size 申请的节点大小
 * \param _nobjs [in, out] 最终申请到的节点个数
 */
template <bool __threads, int __inst>
char* __default_alloc_template<__threads, __inst>::_S_chunk_alloc(size_t _size, int& _nobjs) {
  char* _result = 0;
  size_t _total_size = _size * _nobjs;
  size_t _bytes_left = _S_free_end - _S_free_start;
  if (_bytes_left >= _total_size) {  // 内存池剩余空间完全满足申请
    _result = _S_free_start;
    _S_free_start += _total_size;
    return _result;
  } else if (_bytes_left >= _size) {  // 内存池剩余空间不能满足申请，但是能够提供一个以上的 _size
    // 计算剩余空间能够满足分配的节点数
    _nobjs = (int)(_bytes_left / _size);
    _total_size = _nobjs * _size;
    _result = _S_free_start;
    _S_free_start += _total_size;
    return _result;
  } else {  // 内存池剩余空间无法提供至少一个节点的空间
    // 向堆内存申请空间，申请的空间大小是所要求的2倍还要多一点（申请空间一定是 8 的倍数）
    size_t _bytes_to_get = 2 * _total_size + _S_round_up(_S_heap_size >> 4);
    // 将内存池的剩余空间分给合适的空闲链表
    if (_bytes_left > 0) {
      __Obj* volatile* _my_free_list = _S_free_list + _S_free_list_index(_bytes_left);
      // 将内存池中剩余空间分配给 _my_free_list
      ((__Obj*)_S_free_start)->_next = *_my_free_list;
      *_my_free_list = (__Obj*)_S_free_start;
    }
    // 此时内存池中没有多余空间，那么需要补充内存池
    _S_free_start = (char*)malloc(_bytes_to_get);
    if (0 == _S_free_start) {
      // heap 空间不足，malloc() 失败
      // 由于在 free list 中可能存在可以使用的内存，所以我们会收集那些管理超过 _size 大小内存的 free list 中的内存
      //（不会收集比 _size 小的 free list 中的空间），将其扩充到空闲内存池中
      size_t _i = 0;
      __Obj* volatile* _my_free_list = 0;
      for (_i = _size; _i < _S_max_bytes; _i += _S_align_bytes) {
        _my_free_list = _S_free_list + _S_free_list_index(_i);
        __Obj* _first_node = *_my_free_list;
        if (0 != _first_node) {  // 当前 free list 中尚有未使用的内存
          // 将 _first_node 中内存添加到内存池
          *_my_free_list = _first_node->_next;
          _S_free_start = _first_node->_data;
          _S_free_end = _S_free_start + _i;
          return _S_chunk_alloc(_size, _nobjs);
        }
      }
      _S_free_end = 0;
      // 如果依旧无法扩充内存，则最后使用一级内存配置器
      _S_free_start = (char*)base::malloc_alloc::allocate(_bytes_to_get);
    }
    _S_heap_size += _bytes_to_get;
    _S_free_end = _S_free_start + _bytes_to_get;
    return _S_chunk_alloc(_size, _nobjs);  // 扩充内存池后，递归调用
  }
}

/*!
 * \brief 如果管理大小为 _N 的 free list 中没有可用的节点，则该函数会为该 free list 填充多个新的节点
 * \param _N 指定 free list 管理的节点大小
 */
template <bool __threads, int __inst>
void* __default_alloc_template<__threads, __inst>::_S_refill(size_t _N) {
  int _nobjs = 20;  // 默认一次分配 20 个 _N 大小的块
  char* _chunk = _S_chunk_alloc(_N, _nobjs);
  __Obj* volatile* _my_free_list = 0;
  __Obj* _next_node = 0;
  __Obj* _current_node = 0;
  __Obj* _result = 0;
  int _i = 1;

  // 如果只获得一个节点，那么这个节点就直接分给调用者，自由链表中不会增加新节点
  if (1 == _nobjs) {
    return _chunk;
  }
  // 获取到多个节点，第 0 个节点的内存分配给调用者，其余节点需要被添加到自由链表中
  _result = (__Obj*)_chunk;
  _my_free_list = _S_free_list + _S_free_list_index(_N);
  // 重新设置自由链表头节点（之前 _my_free_list 链表是空的）
  *_my_free_list = (__Obj*)(_chunk + _N);
  _next_node = (__Obj*)(_chunk + _N);
  for (_i = 1; _i < _nobjs; ++_i) {
    _current_node = _next_node;
    _next_node = (__Obj*)(_next_node->_data + _N);
    if (_nobjs - 1 == _i) {
      _current_node->_next = 0;
    } else {
      _current_node->_next = _next_node;
    }
  }
  return _result;
}

typedef __default_alloc_template<__NODE_ALLOCATOR_THREADS, 0> default_alloc;
typedef __default_alloc_template<false, 0> single_thread_alloc;

#ifdef __WIN32_THREADS
// 如果使用 windows 多线程库，这里创建一个实例，仅仅用于初始化内存配置器中的 _S_node_allocator_lock
static base::default_alloc _S_node_allocator_dummy_instance
#endif  // __WIN32_THREADS

}  // namespace base

#endif  // __X_ALLOC_H__
