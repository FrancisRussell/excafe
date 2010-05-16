#ifndef SIMPLE_CFD_NUMERIC_VALUE_MAP_HPP
#define SIMPLE_CFD_NUMERIC_VALUE_MAP_HPP

#include <map>

namespace cfd
{

namespace detail
{

template<typename K, typename V>
class ValueMap
{
public:
  typedef K key_type;
  typedef V value_type;
  typedef std::map<key_type, value_type> internal_map_t;

private:
  internal_map_t map;

public:
  void bind(const key_type& var, const double s)
  {
    map[var] = s;
  }

  const internal_map_t& getReference() const
  {
    return map;
  }
};

}

}

#endif
