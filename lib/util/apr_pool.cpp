#include <apr_pools.h>
#include <excafe/util/apr_pool.hpp>
#include <excafe/util/apr_manager.hpp>

namespace excafe
{

namespace util
{

APRPool::APRPool()
{
  if (!APRManager::instance().isInitialised())
    APRManager::instance().init();

  APRManager::checkSuccess(apr_pool_create(&pool, NULL));
}

APRPool::APRPool(const APRPool& p)
{
  APRManager::checkSuccess(apr_pool_create(&pool, p.pool));
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
