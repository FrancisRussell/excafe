#include <string>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/path.hpp>
#include <simple_cfd/config.h>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/codegen/dynamic_cxx.hpp>

namespace cfd
{

namespace codegen
{

long DynamicCXX::nextID = 0;

fs::path DynamicCXX::getTemp()
{
  // P_tmpdir is (possibly) defined in stdio.h.
  #ifdef P_tmpdir
  return P_tmpdir;
  #else
  return "/tmp";
  #endif
}

void DynamicCXX::writeSource(const fs::path& path) const
{
  fs::ofstream file(path);
  file << code;
  file.close();
}

void DynamicCXX::compileCXX(const fs::path& source, const fs::path& object) const
{
  const std::string sourceString = source.string();
  const std::string objectString = object.string();
  const char* args[] = 
    {EXCAFE_CXX_COMPILER, "-O2", "-shared", "-fPIC", sourceString.c_str(), "-o", objectString.c_str(), NULL};

  const pid_t pid = fork();
  if (pid == 0)
  {
    execvp(args[0], const_cast<char**>(args));
    std::cerr << "exec() of " << args[0] << " failed." << std::endl;
    exit(EXIT_FAILURE);
  }
  else if (pid < 0)
  {
    CFD_EXCEPTION("Call to fork() failed.");
  }
  else
  {
    const int options = 0;
    int status;

    waitpid(pid, &status, options);
    if (!WIFEXITED(status))
    {
      CFD_EXCEPTION("Child process exited abnormally.");
    }
    else if (WEXITSTATUS(status) != EXIT_SUCCESS)
    {
      CFD_EXCEPTION("Child process returned a non-success exit status.");
    }
  }
}

fs::path DynamicCXX::getSourcePath() const
{
  std::ostringstream nameStream;
  nameStream << "excafe_generated_cxx_" << id << ".cpp";
  const fs::path path = getTemp() / nameStream.str();
  return path;
}

fs::path DynamicCXX::getObjectPath() const
{
  std::ostringstream nameStream;
  nameStream << "excafe_generated_lib_" << id << ".so";
  const fs::path path = getTemp() / nameStream.str();
  return path;
}

void DynamicCXX::compile()
{
  const fs::path sourcePath = getSourcePath();
  const fs::path objectPath = getObjectPath();
  writeSource(sourcePath);
  compileCXX(sourcePath, objectPath);
}

}

}
