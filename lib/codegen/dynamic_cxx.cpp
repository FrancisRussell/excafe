#include <string>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_pointer.hpp>
#include <apr_file_io.h>
#include <apr_dso.h>
#include <apr_thread_proc.h>
#include <simple_cfd/config.h>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/util/apr_pool.hpp>
#include <simple_cfd/util/apr_manager.hpp>
#include <simple_cfd/codegen/dynamic_cxx.hpp>

namespace cfd
{

namespace codegen
{

using util::APRManager;
using util::APRPool;

long DynamicCXX::nextID = 0;

DynamicCXX::DynamicCXX(const std::string& _code) : 
  id(nextID++), sourcePath(constructSourcePath()), objectPath(constructObjectPath()),
  compiled(false), code(_code), dsoHandle(NULL)
{
}

std::string DynamicCXX::getTemp() const
{
  APRPool pool;
  const char* temp;
  APRManager::checkSuccess(apr_temp_dir_get(&temp, pool));
  return temp;
}

void DynamicCXX::writeSource(const std::string& path)
{
  std::ofstream file(path.c_str());
  file << code;
  file.close();
}

void DynamicCXX::compile()
{
  compiled = true;
  writeSource(sourcePath);

  const char* args[] = 
    {"c++", "-O2", "-shared", "-fPIC", sourcePath.c_str(), "-o", objectPath.c_str(), NULL};

  APRPool pool;

  apr_procattr_t* attr;
  APRManager::checkSuccess(apr_procattr_create(&attr, pool));
  APRManager::checkSuccess(apr_procattr_cmdtype_set(attr, APR_PROGRAM_PATH));

  apr_proc_t process;
  APRManager::checkSuccess(apr_proc_create(&process, args[0], args, NULL, attr, pool)); 

  int status;
  apr_exit_why_e why;
  apr_proc_wait(&process, &status, &why, APR_WAIT);

  if (why != APR_PROC_EXIT)
  {
    CFD_EXCEPTION("Compiler exited abnormally.");
  }
  else if (status != EXIT_SUCCESS)
  {
    CFD_EXCEPTION("Compiler signalled an error.");
  }
}

void DynamicCXX::loadDSO()
{
  APRManager::checkSuccess(apr_dso_load(&dsoHandle, objectPath.c_str(), classPool));
}

void DynamicCXX::compileAndLoad()
{
  compile();
  loadDSO();
}

apr_dso_handle_sym_t DynamicCXX::getSymbol(const std::string& name)
{
  apr_dso_handle_sym_t symbol;
  APRManager::checkSuccess(apr_dso_sym(&symbol, dsoHandle, name.c_str()));
  return symbol;
}

template<typename T> 
static T convertHandle(const apr_dso_handle_sym_t symbol)
{
  // This function peforms the conversion between function pointers and
  // data pointers needed for using dlsym-style APIs. Rather than
  // casting directly, we dereference a "type-punned" pointer.

  BOOST_STATIC_ASSERT(boost::is_pointer<T>::value);
  BOOST_STATIC_ASSERT(sizeof(T) == sizeof(apr_dso_handle_sym_t));

  const T converted = *reinterpret_cast<T const*>(&symbol);
  return converted;
}

template<> 
apr_dso_handle_sym_t convertHandle(const apr_dso_handle_sym_t symbol)
{
  return symbol;
}

void* DynamicCXX::getData(const std::string& name)
{
  return convertHandle<void*>(getSymbol(name));
}

DynamicCXX::function_t DynamicCXX::getFunction(const std::string& name) 
{
  return convertHandle<function_t>(getSymbol(name));
}

std::string DynamicCXX::mergePath(const std::string& root, const std::string& additional) const
{
  APRPool pool;
  const apr_int32_t flags = 0;
  char* path;
  const apr_status_t result = apr_filepath_merge(&path, root.c_str(), additional.c_str(), flags, pool);
  APRManager::checkSuccess(result);
  return path;
}

std::string DynamicCXX::constructSourcePath() const
{
  std::ostringstream nameStream;
  nameStream << "excafe_generated_cxx_" << id << ".cpp";
  return mergePath(getTemp(), nameStream.str());
}

std::string DynamicCXX::constructObjectPath() const
{
  std::ostringstream nameStream;
  nameStream << "excafe_generated_lib_" << id << ".so";
  return mergePath(getTemp(), nameStream.str());
}

DynamicCXX::~DynamicCXX()
{
  if (dsoHandle != NULL)
    apr_dso_unload(dsoHandle);

  if (compiled)
  {
    remove(sourcePath.c_str());
    remove(objectPath.c_str());
  }
}

}

}
