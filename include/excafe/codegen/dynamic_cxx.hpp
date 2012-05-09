#ifndef EXCAFE_CODEGEN_DYNAMIC_CXX_HPP
#define EXCAFE_CODEGEN_DYNAMIC_CXX_HPP

#include <string>
#include <boost/utility.hpp>
#include <apr_dso.h>
#include <excafe/util/apr_pool.hpp>

namespace excafe
{

namespace codegen
{

class DynamicCXX : public boost::noncopyable
{
private:
  static long nextID;

  const long id;
  const std::string sourcePath;
  const std::string objectPath;

  util::APRPool classPool;
  bool compiled;
  std::string code;
  apr_dso_handle_t* dsoHandle;

  std::string mergePath(const std::string& root, const std::string& additional) const;
  std::string getTemp() const;
  apr_dso_handle_sym_t getSymbol(const std::string& name);
  std::string constructSourcePath() const;
  std::string constructObjectPath() const;
  void writeSource(const std::string& path);
  void compile();
  void loadDSO();

public:
  typedef void (*function_t)();

  DynamicCXX(const std::string& code);
  void compileAndLoad();
  void* getData(const std::string& name);
  function_t getFunction(const std::string& name);
  ~DynamicCXX();
};

}

}

#endif
