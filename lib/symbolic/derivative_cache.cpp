#include <simple_cfd/symbolic/derivative_cache.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/expr.hpp>

namespace cfd
{

namespace symbolic
{
  
DerivativeCache::DerivativeCache() : cache(new cache_t())
{
}

Expr DerivativeCache::derivative(const Expr& e, const Symbol& s)
{
  const cache_t::key_type key(e, s);
  const cache_t::iterator iter = cache->find(key);

  if (iter != cache->end())
  {
    return iter->second;
  }
  else
  {
    const Expr d = e.derivative(s, *this);
    cache->insert(cache_t::value_type(key, d));
    return d;
  }
}

}

}
