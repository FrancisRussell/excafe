#ifndef SIMPLE_CFD_UTIL_EXCEPTION_TRACE_HPP
#define SIMPLE_CFD_UTIL_EXCEPTION_TRACE_HPP

#ifndef __LINUX__
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

class ExceptionTrace
{
private:
  static const std::size_t maxSize = 25;
  int btLength;
  char** symbols;

public:
  ExceptionTrace() : symbols(NULL)
  {
#ifdef __linux__
    void* bt[maxSize];
    btLength = backtrace(bt, maxSize);
    symbols = backtrace_symbols(bt, btLength);

#endif
  }

  const std::string getTrace() const
  {
    if (symbols == NULL)
    {
      return "Stack trace unavailable for this platform.\n";
    }
    else
    {
      std::stringstream traceStream;

      for(int i=0; i<btLength; ++i)
        traceStream << symbols[i] << std::endl;

      return traceStream.str();
    }
  }

  ~ExceptionTrace()
  {
    free(symbols);
  }
};

}

}
#endif
