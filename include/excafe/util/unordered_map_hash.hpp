#ifndef EXCAFE_UTIL_UNORDERED_MAP_HASH_HPP
#define EXCAFE_UTIL_UNORDERED_MAP_HASH_HPP

#include <cstddef>
#include <functional>
#include <excafe/util/hash.hpp>

namespace excafe
{

namespace util
{

template<typename Map>
struct UnorderedMapHash : public std::unary_function<Map, std::size_t>
{
  std::size_t operator()(const Map& map) const
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
};

}

}

#endif
