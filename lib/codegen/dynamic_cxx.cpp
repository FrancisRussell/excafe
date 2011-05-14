#include <string>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <apr_file_io.h>
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
  id(nextID++), code(_code)
{
  APRManager::checkSuccess(apr_procattr_create(&attr, classPool));
  APRManager::checkSuccess(apr_procattr_cmdtype_set(attr, APR_PROGRAM_PATH));
}

std::string DynamicCXX::getTemp() const
{
  APRPool pool;
  const char* temp;
  APRManager::checkSuccess(apr_temp_dir_get(&temp, pool));
  return temp;
}

void DynamicCXX::writeSource(const std::string& path) const
{
  std::ofstream file(path.c_str());
  file << code;
  file.close();
}

void DynamicCXX::compileCXX(const std::string& source, const std::string& object) const
{
  APRPool pool;
  const char* args[] = 
    {"c++", "-O2", "-shared", "-fPIC", source.c_str(), "-o", object.c_str(), NULL};

  apr_proc_t process;
  APRManager::checkSuccess(apr_proc_create(&process, args[0], args, NULL, attr, pool)); 

  int status;
  apr_exit_why_e why;
  apr_proc_wait(&process, &status, &why, APR_WAIT);

  if (why != APR_PROC_EXIT)
  {
    CFD_EXCEPTION("Child process exited abnormally.");
  }
  else if (status != EXIT_SUCCESS)
  {
    CFD_EXCEPTION("Child process returned a non-success exit status.");
  }
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

std::string DynamicCXX::getSourcePath() const
{
  std::ostringstream nameStream;
  nameStream << "excafe_generated_cxx_" << id << ".cpp";
  return mergePath(getTemp(), nameStream.str());
}

std::string DynamicCXX::getObjectPath() const
{
  std::ostringstream nameStream;
  nameStream << "excafe_generated_lib_" << id << ".so";
  return mergePath(getTemp(), nameStream.str());
}

void DynamicCXX::compile()
{
  const std::string sourcePath = getSourcePath();
  const std::string objectPath = getObjectPath();
  writeSource(sourcePath);
  compileCXX(sourcePath, objectPath);
}

}

}
