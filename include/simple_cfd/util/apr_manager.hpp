#ifndef SIMPLE_CFD_UTIL_APR_MANAGER_HPP
#define SIMPLE_CFD_UTIL_APR_MANAGER_HPP

#include <string>
#include <apr_errno.h>
#include "simple_cfd_fwd.hpp"

namespace cfd
{

namespace util
{

class APRManager
{
private:
  friend class util::Singleton<APRManager>;
  bool initialised;

  APRManager();
  APRManager(const APRManager&);
  APRManager& operator=(const APRManager&);

public:
  static APRManager& instance();
  static std::string getDescription(apr_status_t status);
  static void checkSuccess(apr_status_t status);

  void init();
  bool isInitialised() const;
  ~APRManager();
};

}

}

#endif
