#ifndef SIMPLE_CFD_MP_INTEGER_HPP
#define SIMPLE_CFD_MP_INTEGER_HPP

#include "mp_fwd.hpp"
#include <gmp.h>
#include <iosfwd>
#include <climits>
#include <vector>
#include <boost/operators.hpp>
#include <boost/shared_array.hpp>
#include <simple_cfd/util/hybrid_array.hpp>
#include <simple_cfd/util/hash.hpp>
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace mp
{

class Integer : boost::totally_ordered<Integer>,
                boost::totally_ordered<Integer, int>,
                boost::totally_ordered<Integer, long>,
                boost::integer_arithmetic<Integer>,
                boost::shiftable<Integer, unsigned>
{
private:
  friend class Rational;
  friend class Float;

  static const int STACK_LIMBS = 4;
  struct width_tag {};

  int size;
  int allocated;
  boost::shared_array<mp_limb_t> data;

  int width() const
  {
    return size < 0 ? -size : size;
  }

  bool isNegative() const
  {
    return size < 0;
  }

  const mp_limb_t* limbs() const
  {
    return data.get();
  }

  mp_limb_t* limbs()
  {
    return data.get();
  }
  
  static int negate(const bool b, const int val)
  {
    return (b ? -val : val);
  }

  static std::size_t numLimbs(const int bits);
  void reallocUnique(const int newAllocated);
  int computeWidth(const int maxWidth) const;
  void performAddition(const Integer* u, const Integer* v);

  template<typename T>
  void initialise(const T& i)
  {
    reallocUnique(numLimbs(sizeof(T)*CHAR_BIT));

    const bool negative = i<0;
    T magnitude = (negative ? -i : i);
    int loc = 0;
    while(magnitude>0)
    {
      data[loc] = magnitude & GMP_NUMB_MASK;
      magnitude /= GMP_NUMB_MAX;
      ++loc;
    }

    size = negate(negative, loc);
  }

  template<typename T>
  bool equal(const T& i) const
  {
    if (size == 0 || i == 0)
    {
      // Both zero
      return size == 0 && i == 0;
    }
    else if ((i < 0) != isNegative())
    {
      // Differing signs
      return false;
    }

    const T absI = std::abs(i);
    if (absI < GMP_NUMB_MAX)
    {
      // i is too small to be equal
      if (width() > 1)
        return false;

      // Both values fit into a single limb
      const mp_limb_t limb = static_cast<mp_limb_t>(absI);
      return (data[0] & GMP_NUMB_MASK) == limb;
    }

    return *this == Integer(i);
  }

  template<typename T>
  bool lessThan(const T& i) const
  {
    if (size == 0)
    {
      // Any positive i is larger
      return i>0;
    }
    else if ((i < 0) != isNegative())
    {
      // Differing signs
      return isNegative();
    }

    const T absI = std::abs(i);
    if (absI < GMP_NUMB_MAX)
    {
      // i has a smaller magnitude
      if (width() > 1)
        return isNegative();

      // Both values fit into a single limb
      mp_limb_t left = data[0];
      mp_limb_t right = static_cast<mp_limb_t>(absI);

      if (isNegative()) 
        std::swap(left, right);

      return left < right;
    }

    return *this < Integer(i);
  }

  void mpzInit(mpz_t& mpz) const;
  Integer(const int size, const int allocated, const boost::shared_array<mp_limb_t>& data);

public:
  static Integer gcd(const Integer& x, const Integer& y);
  static Integer lcm(const Integer& x, const Integer& y);

  Integer();
  Integer(const char* str);
  Integer(const Integer& i);
  Integer(const int i);
  Integer(const long i);
  Integer(const unsigned int i);
  Integer(const unsigned long i);

  bool operator==(const Integer& i) const;
  bool operator==(const int i) const;
  bool operator==(const long i) const;
  bool operator<(const Integer& i) const;
  bool operator<(const int i) const;
  bool operator<(const long i) const;
  Integer& operator=(const Integer& i);
  Integer& operator+=(const Integer& i);
  Integer& operator-=(const Integer& i);
  Integer& operator++();
  Integer& operator--();
  Integer& operator*=(const Integer& i);
  Integer& operator/=(const Integer& dividend);
  Integer& operator%=(const Integer& dividend);
  Integer& operator>>=(unsigned i);
  Integer& operator<<=(unsigned i);
  Integer operator-() const;

  Integer abs() const;
  Integer isqrt() const;
  Integer pow(const long exponent) const;
  std::size_t hash() const;
  void write(std::ostream& out) const;
  void swap(Integer& i);
  int toInt() const;
  long toLong() const;
  float toFloat() const;
  double toDouble() const;
};

std::ostream& operator<<(std::ostream& o, const Integer& i);
std::size_t hash_value(const Integer& i);
Integer abs(const Integer& i);
Integer pow(const Integer& i, const int exponent);
Integer pow(const Integer& i, const long exponent);
Integer isqrt(const Integer& i);

}

}

#endif
