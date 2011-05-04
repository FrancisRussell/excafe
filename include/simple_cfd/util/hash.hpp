#ifndef SIMPLE_CFD_UTIL_HASH_HPP
#define SIMPLE_CFD_UTIL_HASH_HPP

#include <cstddef>
#include <climits>
#include <boost/functional/hash.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_unsigned.hpp>
#include <cln/integer.h>

/* Based on the FNV-1a mixing function[0]. Performing the XOR before the
   multiply results in superior avalanche behaviour [1].

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

template<typename T>
static void hash_accum_integral(std::size_t& seed, const T value)
{
  BOOST_STATIC_ASSERT(boost::is_integral<T>::value && boost::is_unsigned<T>::value);

  for(std::size_t i=0; i<sizeof(T); ++i)
  {
    const std::size_t shift = i*CHAR_BIT;
    seed ^= (value >> shift) & ((1 << CHAR_BIT) - 1);
    seed *= detail::fnv_traits<sizeof(std::size_t)*CHAR_BIT>::prime;
  }
}

}

template<typename T>
static void hash_accum(std::size_t& seed, const T& value)
{
  const std::size_t valueHash = boost::hash<T>()(value);
  detail::hash_accum_integral(seed, valueHash);
}

template<>
void hash_accum(std::size_t& seed, const cln::cl_I& value)
{
  static const std::size_t uLongWidth = CHAR_BIT*sizeof(unsigned long);

  // We add one for the sign bit.
  const std::size_t width = cln::integer_length(value) + 1;
  for(std::size_t start=0; start<width; start+=uLongWidth)
  {
    const cln::cl_byte bits(uLongWidth, start);
    const unsigned long next = cln::cl_I_to_ulong(cln::ldb(value, bits));
    detail::hash_accum_integral(seed, next);
  }
}

}

}

#endif
