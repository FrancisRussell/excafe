#ifndef SIMPLE_CFD_UTIL_APR_MANAGER_HPP
#define SIMPLE_CFD_UTIL_APR_MANAGER_HPP

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
  void init();
  bool isInitialised() const;
  ~APRManager();
};

}

}

#endif
