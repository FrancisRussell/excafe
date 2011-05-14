#ifndef SIMPLE_CFD_CODEGEN_DYNAMIC_CXX_HPP
#define SIMPLE_CFD_CODEGEN_DYNAMIC_CXX_HPP

#include <string>
#include <boost/utility.hpp>
#include <simple_cfd/util/apr_pool.hpp>
#include <apr_thread_proc.h>

namespace cfd
{

namespace codegen
{

class DynamicCXX : public boost::noncopyable
{
private:
  static long nextID;

  util::APRPool classPool;
  apr_procattr_t* attr;
  long id;
  std::string code;

  std::string mergePath(const std::string& root, const std::string& additional) const;
  
private:
  std::string getTemp() const;
  void writeSource(const std::string& path) const;
  void compileCXX(const std::string& source, const std::string& object) const;

public:
  DynamicCXX(const std::string& _code);
  std::string getSourcePath() const;
  std::string getObjectPath() const;
  void compile();
};

}

}

#endif
