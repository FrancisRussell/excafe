#ifndef EXCAFE_UTIL_HASH_UNORDERED_MAP_HPP
#define EXCAFE_UTIL_HASH_UNORDERED_MAP_HPP

#include <cstddef>
#include <excafe/util/hash.hpp>

namespace excafe
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
    excafe::util::hash_accum(termHash, iter->first);
    excafe::util::hash_accum(termHash, iter->second);
    result ^= termHash;
  }

  excafe::util::hash_accum(result, map.size());
  return result;
}

}

}

#endif
