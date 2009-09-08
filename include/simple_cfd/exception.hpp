#ifndef SIMPLE_CFD_EXCEPTION_HPP
#define SIMPLE_CFD_EXCEPTION_HPP

#include <exception>
#include <sstream>
#include "util/exception_trace.hpp"

#define CFD_EXCEPTION(x) throw CFDException(x, __FILE__, __LINE__)

namespace cfd
{

class CFDException : public std::exception
{
private:
  const util::ExceptionTrace trace;
  const std::string description;
  const std::string file;
  const int line;
  mutable std::string whatString;

public:
  CFDException(const std::string _description, const char* const _file, const int _line) : 
    description(_description), file(_file), line(_line)
  {
  }

  virtual const char* what() const throw()
  {
    std::stringstream whatStream;
    whatStream << file << ":" << line << ": " << description << std::endl;
    whatStream << "Stack trace: " << std::endl << trace.getTrace();
    whatString = whatStream.str();

    return whatString.c_str();
  }

  virtual ~CFDException() throw()
  {
  }
};

}

#endif
