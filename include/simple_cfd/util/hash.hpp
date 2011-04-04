#ifndef SIMPLE_CFD_UTIL_HASH_HPP
#define SIMPLE_CFD_UTIL_HASH_HPP

#include <cstddef>
#include <climits>
#include <boost/functional/hash.hpp>
#include <boost/static_assert.hpp>

// See http://www.isthe.com/chongo/tech/comp/fnv/index.html for details
// on the hashing function used.

namespace cfd
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
    seed *= detail::fnv_traits<sizeof(std::size_t)*CHAR_BIT>::prime;
    seed ^= valueHash & ((1 << CHAR_BIT) - 1);
    valueHash >>= CHAR_BIT;
  }
}


}

#endif
