#ifndef SIMPLE_CFD_CODEGEN_DYNAMIC_CXX_HPP
#define SIMPLE_CFD_CODEGEN_DYNAMIC_CXX_HPP

#include <string>
#include <boost/filesystem/path.hpp>

namespace cfd
{

namespace codegen
{

namespace fs = boost::filesystem; 

class DynamicCXX
{
private:
  static long nextID;

  long id;
  std::string code;
  
private:
  static fs::path getTemp();
  void writeSource(const fs::path& path) const;
  void compileCXX(const fs::path& source, const fs::path& object) const;

public:
  DynamicCXX(const std::string& _code) : 
    id(nextID++), code(_code)
  {
  }

  fs::path getSourcePath() const;
  fs::path getObjectPath() const;
  void compile();
};

}

}

#endif
