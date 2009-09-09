#ifndef SIMPLE_CFD_UTIL_EXCEPTION_TRACE_HPP
#define SIMPLE_CFD_UTIL_EXCEPTION_TRACE_HPP

#if defined(__linux__)
extern "C" {
#include <execinfo.h>
}
#endif

#include <string>
#include <sstream>

namespace cfd
{

namespace util
{

#if defined(__linux__)

class LinuxExceptionTrace
{
private:
  static const std::size_t maxSize = 25;
  void* bt[maxSize];
  int btLength;

public:
  LinuxExceptionTrace()
  {
    btLength = backtrace(bt, maxSize);
  }

  const std::string getTrace() const
  {
    char** const symbols = backtrace_symbols(bt, btLength);
    std::stringstream traceStream;

    for(int i=0; i<btLength; ++i)
      traceStream << symbols[i] << std::endl;

    free(symbols);
    return traceStream.str();
  }
};
typedef LinuxExceptionTrace ExceptionTrace;

#else

class GenericExceptionTrace
{
public:
  GenericExceptionTrace()
  {
  }

  const std::string getTrace() const
  {
    return "Stack trace unavailable for this platform.\n";
  }
};
typedef GenericExceptionTrace ExceptionTrace;

#endif

}

}
#endif
