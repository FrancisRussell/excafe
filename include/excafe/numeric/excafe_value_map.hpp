#ifndef EXCAFE_NUMERIC_EXCAFE_VALUE_MAP_HPP
#define EXCAFE_NUMERIC_EXCAFE_VALUE_MAP_HPP

#include <excafe/symbolic/expr.hpp>
#include "excafe_mapper.hpp"

namespace excafe
{

namespace detail
{

template<typename K, typename V>
class ExcafeValueMap
{
public:
  typedef K key_type;
  typedef V value_type;
  typedef symbolic::Expr::subst_map internal_map_t;

private:
  internal_map_t map;

public:
  void bind(const key_type& var, const double s)
  {
    detail::ExcafeMapper<key_type>& mapper(detail::ExcafeMapper<key_type>::instance());
    map[mapper.getSymbol(var)] = s;
  }

  void bind(const key_type& var, const value_type& value)
  {
    detail::ExcafeMapper<key_type>& mapper(detail::ExcafeMapper<key_type>::instance());
    map[mapper.getSymbol(var)] = value.expr;
  }

  const internal_map_t& getReference() const
  {
    return map;
  }
};

}

}

#endif
