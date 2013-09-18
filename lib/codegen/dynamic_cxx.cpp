#include <string>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <sstream>
#include <fstream>
#include <boost/static_assert.hpp>
#include <boost/integer.hpp>
#include <boost/type_traits/is_pointer.hpp>
#include <apr_file_io.h>
#include <apr_dso.h>
#include <apr_thread_proc.h>
#include <apr_user.h>
#include <excafe/config.h>
#include <excafe/exception.hpp>
#include <excafe/util/apr_pool.hpp>
#include <excafe/util/apr_manager.hpp>
#include <excafe/codegen/dynamic_cxx.hpp>

namespace excafe
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
  if (dsoHandle == NULL)
    CFD_EXCEPTION("Cannot get symbol because DSO is not loaded!");

  apr_dso_handle_sym_t symbol;
  APRManager::checkSuccess(apr_dso_sym(&symbol, dsoHandle, name.c_str()));
  return symbol;
}

template<typename T>
static T convertHandle(const apr_dso_handle_sym_t symbol)
{
  // This function peforms the conversion between function pointers and
  // data pointers needed for using dlsym-style APIs. We cast via an
  // integral type of the same size.

  BOOST_STATIC_ASSERT(boost::is_pointer<T>::value);
  BOOST_STATIC_ASSERT(boost::is_pointer<apr_dso_handle_sym_t>::value);
  BOOST_STATIC_ASSERT(sizeof(T) == sizeof(apr_dso_handle_sym_t));

  typedef typename boost::int_t<CHAR_BIT*sizeof(T)>::exact integral_t;
  const integral_t integral = reinterpret_cast<integral_t>(symbol);
  const T converted = reinterpret_cast<T>(integral);
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

std::string DynamicCXX::getNamePrefix() const
{
  std::ostringstream prefixStream;
  prefixStream << "excafe_generated";

#ifdef APR_HAS_USER
  APRPool pool;
  apr_uid_t uid;
  apr_gid_t gid;

  const apr_status_t result = apr_uid_current(&uid, &gid, pool);

  if (result == APR_SUCCESS)
    prefixStream << "_" << uid;
#endif

  return prefixStream.str();
}

std::string DynamicCXX::constructSourcePath() const
{
  std::ostringstream nameStream;
  nameStream << getNamePrefix() << "_cxx_" << id << ".cpp";
  return mergePath(getTemp(), nameStream.str());
}

std::string DynamicCXX::constructObjectPath() const
{
  std::ostringstream nameStream;
  nameStream << getNamePrefix() << "_lib_" << id << ".so";
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
