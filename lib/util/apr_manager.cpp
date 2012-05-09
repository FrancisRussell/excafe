#include <string>
#include <vector>
#include <apr_general.h>
#include <excafe/util/singleton.hpp>
#include <excafe/util/apr_manager.hpp>
#include <excafe/exception.hpp>

namespace excafe
{

namespace util
{

APRManager::APRManager() : initialised(false)
{
}

APRManager& APRManager::instance()
{
  return util::Singleton<APRManager>::getInstance();
}

void APRManager::init()
{
  if (!initialised)
  {
    checkSuccess(apr_initialize());
    initialised = true;
  }
  else
  {
    CFD_EXCEPTION("Attempted to initialise APRManager multiple times.");
  }
}

bool APRManager::isInitialised() const
{
  return initialised;
}

void APRManager::checkSuccess(const apr_status_t status)
{
  if (status != APR_SUCCESS)
    CFD_EXCEPTION(getDescription(status));
}

std::string APRManager::getDescription(const apr_status_t status)
{
  std::vector<char> buffer(256);
  apr_strerror(status, &buffer[0], buffer.size());
  return std::string(&buffer[0]);
}
  
APRManager::~APRManager()
{
  apr_terminate();
}

}

}
