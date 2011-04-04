#ifndef SIMPLE_CFD_UTIL_HASH_UNORDERED_MAP_HPP
#define SIMPLE_CFD_UTIL_HASH_UNORDERED_MAP_HPP

#include <cstddef>
#include <simple_cfd/util/hash.hpp>

namespace cfd
{

namespace util
{

template<typename Map>
std::size_t hash_unordered_map(const Map& map)
{
  std::size_t result = 0x2cb6bfcd;
  for(typename Map::const_iterator iter = map.begin();
      iter != map.end();
      ++iter)
  {   
    std::size_t termHash = 0x2b20b3a2;
    cfd::hash_accum(termHash, *iter);
    result ^= termHash;
  }

  cfd::hash_accum(result, map.size());
  return result;
}

}

}

#endif
