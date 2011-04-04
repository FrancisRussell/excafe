#ifndef SIMPLE_CFD_UTIL_HASH_HPP
#define SIMPLE_CFD_UTIL_HASH_HPP

#include <cstddef>
#include <climits>
#include <boost/functional/hash.hpp>
#include <boost/static_assert.hpp>

/* Based on the FNV-1a mixing function[0]. Performing the XOR before the
   multiply results in much better avalanche behaviour (apparently)
   [1].

   [0] http://www.isthe.com/chongo/tech/comp/fnv/index.html 
   [1] https://sites.google.com/site/murmurhash/avalanche
*/

namespace cfd
{

namespace util
{

namespace detail
{

template<std::size_t width>
struct fnv_traits
{
};

template<>
struct fnv_traits<32>
{
  static const unsigned long basis = 2166136261ul;
  static const unsigned long prime = 16777619ul;
};

template<>
struct fnv_traits<64>
{
  static const unsigned long long basis = 14695981039346656037ull;
  static const unsigned long long prime = 1099511628211ull;
};

}

template<typename T>
void hash_accum(std::size_t& seed, const T& value)
{
  std::size_t valueHash = boost::hash<T>()(value);
  for(std::size_t i=0; i<sizeof(std::size_t); ++i)
  {
    seed ^= valueHash & ((1 << CHAR_BIT) - 1);
    seed *= detail::fnv_traits<sizeof(std::size_t)*CHAR_BIT>::prime;
    valueHash >>= CHAR_BIT;
  }
}

}

}

#endif
