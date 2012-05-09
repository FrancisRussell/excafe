#ifndef EXCAFE_UTIL_APR_POOL_HPP
#define EXCAFE_UTIL_APR_POOL_HPP

#include <apr_pools.h>
#include <boost/utility.hpp>

namespace excafe
{

namespace util
{

class APRPool : public boost::noncopyable
{
private:
  apr_pool_t* pool;

  APRPool& operator=(const APRPool&);

public:
  APRPool();
  APRPool(const APRPool&);
  operator apr_pool_t*() const;
  ~APRPool();
};

}

}

#endif
