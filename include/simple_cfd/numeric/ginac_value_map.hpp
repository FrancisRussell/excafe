#ifndef SIMPLE_CFD_NUMERIC_GINAC_VALUE_MAP_HPP
#define SIMPLE_CFD_NUMERIC_GINAC_VALUE_MAP_HPP

#include <ginac/basic.h>
#include "ginac_mapper.hpp"

namespace cfd
{

namespace detail
{

template<typename K>
class GinacValueMap
{
public:
  typedef K key_type;
  typedef GiNaC::exmap internal_map_t;

private:
  internal_map_t map;

public:
  void bind(const key_type& var, const double s)
  {
    detail::GinacMapper<key_type>& mapper(detail::GinacMapper<key_type>::instance());
    map[mapper.getGiNaCSymbol(var)] = s;
  }

  const internal_map_t& getReference() const
  {
    return map;
  }
};

}

}

#endif
