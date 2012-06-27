#ifndef EXCAFE_UTIL_ORDERED_MAP_HASH_HPP
#define EXCAFE_UTIL_ORDERED_MAP_HASH_HPP

#include <cstddef>
#include <functional>
#include <excafe/util/hash.hpp>

namespace excafe
{

namespace util
{

template<typename Map>
struct OrderedMapHash : public std::unary_function<Map, std::size_t>
{
  std::size_t operator()(const Map& map) const
  {
    std::size_t result = 0xbf51b047;
    std::size_t size = 0;
    for(typename Map::const_iterator iter = map.begin();
        iter != map.end();
        ++iter)
    { 
      ++size;
      excafe::util::hash_accum(result, 0xa1c6605d);
      excafe::util::hash_accum(result, iter->first);
      excafe::util::hash_accum(result, iter->second);
    }

    excafe::util::hash_accum(result, size);
    return result;
  }
};

}

}

#endif
