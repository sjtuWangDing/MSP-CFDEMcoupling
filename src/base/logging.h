#ifndef __LOGGING_H__
#define __LOGGING_H__

#include <inttypes.h>
#include <time.h>
#include <stdarg.h> // va_list
#include <iostream>
#include <cstdlib>  // abort()
#include <sys/types.h>
#include <unistd.h>
#include <pthread.h>

#define NONE "\033[m"
#define RED "\033[1;31m"
#define GREEN "\033[1;32m"
#define YELLOW "\033[1;33m"

namespace base {

class Logger {
public:
  //! \brief 日志等级 - 调试
  static const int LEVEL_DBG = 0;
  //! \brief 日志等级 - 信息
  static const int LEVEL_INF = 1;
  //! \brief 日志等级 - 警告
  static const int LEVEL_WAR = 2;
  //! \brief 日志等级 - 一般错误
  static const int LEVEL_ERR = 3;
  //! \brief 日志等级 - 致命错误
  static const int LEVEL_CRT = 4;
  /*!
   * \brief 日志等级标签数组
   *        debug   - 调试
   *        info    - 信息
   *        warning - 警告
   *        error   - 一般错误
   *        crt     - 致命错误
   */
  static constexpr const char *LEVELS[] = {
    "debug",
    "info",
    "warning",
    "error",
    "crt"
  };

  //! \brief 获取文件名
  inline static std::string getFileName(const std::string& path) {
    char ch = '/';
#ifdef _WIN32
    ch = '\\';
#endif // _WIN32
    std::string::size_type pos = path.rfind(ch);
    if (std::string::npos == pos) {
      return path;
    } else {
      return path.substr(pos + 1);
    }
  }

  /*!
   * \brief print log according to format
   * \param filePath 源码文件路径
   * \param line 源码行号
   * \param format 格式
   * \param logLevel 日志等级
   */
  inline static void traceDebug(const char* filePath,
                                int line,
                                const char* format,
                                int logLevel,
                                ...) {
//    // 格式化当期系统日期和时间字符串
//    char dataTime[128];
//    time_t now = time((time_t*)nullptr);
//    strftime(dataTime, sizeof(dataTime), "%Y-%m-%d %H:%M:%S", localtime(&now));
//
//    fprintf(stdout, "[%s][%s][pid = %u][tid = %lu][%s : %d]\n",
//            dataTime, LEVELS[logLevel], getpid(), pthread_self(), getFileName(filePath).c_str(), line);
//    switch (logLevel) {
//      case LEVEL_CRT:
//        fprintf(stdout, RED "%s: " NONE, "Critical Error");
//        break;
//      case LEVEL_ERR:
//        fprintf(stdout, RED "%s: " NONE, "Error");
//        break;
//      case LEVEL_WAR:
//        fprintf(stdout, YELLOW "%s: " NONE, "Warming");
//        break;
//      case LEVEL_DBG:
//        fprintf(stdout, GREEN "%s: " NONE, "Debug");
//        break;
//      default:
//        fprintf(stdout, NONE "%s: " NONE, "Info");
//    }
//    va_list args; // 变长参数表
//    va_start(args, format); // 用format以后的参数初始化变长参数表
//    vfprintf(stdout, format, args); // 按format格式打印变长参数表中的内容
//    va_end(args); // 销毁变长参数表
//    fprintf (stdout, "\n\n");
  }
};

#ifndef __TRACE_DEBUG
#define __TRACE_DEBUG(...) base::Logger::traceDebug(__FILE__, __LINE__, __FUNCION__, __VA_ARGS__);
#endif // __TRACE_DEBUG

class LogMessage {
public:
  LogMessage(const char* file, int line): log_stream_(std::cerr) {
    log_stream_ << file << ": " << line << ": ";
  }

  ~LogMessage() {
    log_stream_ << std::endl;
  }

  inline std::ostream& stream() { return log_stream_; }

protected:
  std::ostream& log_stream_;

private:
  LogMessage(const LogMessage&);
  LogMessage& operator= (const LogMessage&);
};

class LogMessageFatal: public LogMessage {
public:
  LogMessageFatal(const char* file, int line): LogMessage(file, line) {}

  ~LogMessageFatal() {
    abort();
  }

private:
  LogMessageFatal(const LogMessageFatal&);
  LogMessageFatal& operator= (const LogMessageFatal&);
};

/*! \brief double type that will be used in default by cfdem */
typedef double d_real_t;

/*! \brief device name CPU */
struct cpu {
  /*! \brief whether this device is CPU or not */
  static const bool kDevCPU = true;
  /*! \brief device flag number, identifies this device */
  static const int kDevMask = 1 << 0;
};

/*! \brief device name GPU */
struct gpu {
  /*! \brief whether this device is CPU or not */
  static const bool kDevCPU = false;
  /*! \brief device flag number, identifies this device */
  static const int kDevMask = 1 << 1;
};

} // namespace base

/*! \brief type that will be used for index */
#ifndef CFDEM_INT64_TENSOR_SIZE
#define CFDEM_INT64_TENSOR_SIZE 0
#endif

/*! \brief whether do padding during allocation */
#ifndef CFDEM_ALLOC_PAD
#define CFDEM_ALLOC_PAD false
#endif

/*! \brief whether use SSE */
#ifndef CFDEM_USE_SSE
#define CFDEM_USE_SSE 1
#endif

#if CFDEM_INT64_TENSOR_SIZE == 1
typedef int64_t index_t;
#else
typedef int32_t index_t;
#endif

#define CHECK(x)                                                                                           \
if (!(x))                                                                                                  \
base::LogMessageFatal(base::Logger::getFileName(__FILE__).c_str(), __LINE__).stream()                      \
  << "Check (" << #x << ") faild! "
#define CHECK_GT(x, y) CHECK((x) > (y))
#define CHECK_LE(x, y) CHECK((x) <= (y))
#define CHECK_GE(x, y) CHECK((x) >= (y))
#define CHECK_EQ(x, y) CHECK((x) == (y))
#define CHECK_NE(x, y) CHECK((x) != (y))
#define CHECK_NOTNULL(x)                                                                                   \
((x) == NULL ? base::LogMessageFatal(base::Logger::getFileName(__FILE__).c_str(), __LINE__).stream()       \
  << "Check not null failed: " #x, (x) : (x))

#endif // __LOGGING_H__
