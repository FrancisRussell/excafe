#include <apr_pools.h>
#include <simple_cfd/util/apr_pool.hpp>
#include <simple_cfd/util/apr_manager.hpp>

namespace cfd
{

namespace util
{

APRPool::APRPool()
{
  if (!APRManager::instance().isInitialised())
    APRManager::instance().init();

  apr_pool_create(&pool, NULL);
}

APRPool::APRPool(const APRPool& p)
{
  apr_pool_create(&pool, p.pool);
}

APRPool::operator apr_pool_t*() const
{
  return pool;
}

APRPool::~APRPool()
{
  apr_pool_destroy(pool);
}

}

}
