#ifndef SIMPLE_CFD_UTIL_APR_POOL_HPP
#define SIMPLE_CFD_UTIL_APR_POOL_HPP

#include <apr_pools.h>
#include <boost/utility.hpp>

namespace cfd
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
