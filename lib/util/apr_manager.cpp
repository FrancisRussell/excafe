#include <apr_general.h>
#include <simple_cfd/util/singleton.hpp>
#include <simple_cfd/util/apr_manager.hpp>
#include <simple_cfd/exception.hpp>

namespace cfd
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
    apr_initialize();
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
  
APRManager::~APRManager()
{
  apr_terminate();
}

}

}
