#ifndef SIMPLE_CFD_MP_INTEGER_HPP
#define SIMPLE_CFD_MP_INTEGER_HPP

#include "mp_fwd.hpp"
#include <gmp.h>
#include <apr_general.h>
#include <iosfwd>
#include <climits>
#include <cstring>
#include <cmath>
#include <stdint.h>
#include <vector>
#include <cassert>
#include <boost/operators.hpp>
#include <boost/mpl/min_max.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/integer_traits.hpp>
#include <boost/type_traits/make_unsigned.hpp>
#include <boost/static_assert.hpp>

#ifdef EXCAFE_VALIDATE_MP
#define EXCAFE_VALIDATE_INTEGER this->validate();
#else
#define EXCAFE_VALIDATE_INTEGER
#endif

namespace excafe
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

  static const int TAG_BITS = 1;
  static const int VALUE_TAG = 1;

  static const Integer ONE;
  static const Integer MINUS_ONE;

  static const mp_limb_t MAX_PACKED_LIMB = boost::mpl::min<
      boost::mpl::integral_c<mp_limb_t, GMP_NUMB_MAX>,
      boost::mpl::integral_c<uintptr_t, (boost::integer_traits<uintptr_t>::const_max >> TAG_BITS)>
    >::type::value;

  struct Header
  {
    int allocated;
    apr_uint32_t count;
  };

  class ConstPacker
  {
  private:
    ConstPacker(const ConstPacker&);
    ConstPacker& operator=(const ConstPacker&);

    const Integer* integer;

  protected:
    mp_limb_t  local[2];
    mp_limb_t* data;

    inline bool isPacked() const
    {
      return integer->data & VALUE_TAG;
    }

    inline int getAllocated() const
    {
      if (isPacked())
        return 2;
      else
        return countLimbs(data);
    }

  public:
    void unpack()
    {
      if (isPacked())
      {
        const uintptr_t value = (integer->data >> TAG_BITS);
        local[0] = value;
        data = local;
      }
      else
      {
        data = reinterpret_cast<mp_limb_t*>(integer->data);
      }
    }

    ConstPacker(const Integer* _integer) : integer(_integer)
    {
      unpack();
    }

    const mp_limb_t* limbs() const
    {
      return data;
    }

    int computeWidth(const int maxWidth) const
    {
      for(int cwidth = maxWidth; cwidth>0; --cwidth)
      {
        if (limbs()[cwidth-1] != 0)
          return cwidth;
      }
      return 0;
    }
  };

  class Packer : public ConstPacker
  {
  private:
    Integer* integer;
    bool complete;

  public:
    Packer(Integer* _integer) : ConstPacker(_integer), integer(_integer), complete(false)
    {
    }

    mp_limb_t* limbs()
    {
      return data;
    }

    void reallocUnique(const int newAllocated)
    {
      const bool packed = isPacked();

      if ((packed && newAllocated > 2)
          || (!packed && (newAllocated > getAllocated() || !isUnique(limbs()))))
      {
        mp_limb_t* const newLimbs = allocateLimbs(newAllocated);
        memcpy(newLimbs, limbs(), std::min(getAllocated(), newAllocated) * sizeof(mp_limb_t));

        if (!packed)
          decrementUseCount(limbs());

        integer->data = reinterpret_cast<uintptr_t>(newLimbs);
        data = newLimbs;
      }
    }

    // The size of the result *must* be set before this is called, otherwise the resulting commit
    // may truncate the result
    void commit()
    {
      assert(!complete);
      complete = true;

      if (isPacked())
      {
        if (integer->width() < 2 && local[0] <= MAX_PACKED_LIMB)
        {
          integer->data = (static_cast<uintptr_t>(local[0]) << TAG_BITS) | VALUE_TAG;
        }
        else
        {
          mp_limb_t* const newLimbs = allocateLimbs(2);
          newLimbs[0] = local[0];
          newLimbs[1] = local[1];
          integer->data = reinterpret_cast<uintptr_t>(newLimbs);
          data = newLimbs;
        }
      }
    }

    void abort()
    {
      complete = true;
    }

    ~Packer()
    {
      assert(complete && "All non-const packers must be explicitly committed or aborted.");
    }
  };

  static const int STACK_LIMBS = 4;

  int size;
  uintptr_t data;

  int width() const
  {
    return size < 0 ? -size : size;
  }

  inline bool isNegative() const
  {
    return size < 0;
  }

  static inline int negate(const bool b, const int val)
  {
    return (b ? -val : val);
  }

  static std::size_t numLimbs(const int bits);
  static mp_limb_t* allocateLimbs(const std::size_t count);
  static std::size_t countLimbs(const mp_limb_t* count);
  static void incrementUseCount(mp_limb_t* limbs);
  static void decrementUseCount(mp_limb_t* limbs);
  static bool isUnique(mp_limb_t* limbs);

  void validate() const;
  void reallocUnique(const int newAllocated);
  void performAddition(const Integer* u, const Integer* v);

  template<typename T>
  void initialise(T i)
  {
    BOOST_STATIC_ASSERT(boost::integer_traits<T>::is_integral);

    Packer packer(this);
    packer.reallocUnique(numLimbs(sizeof(T)*CHAR_BIT));

    // We need to handle the special case where
    // i is the minimum possible signed value.
    const bool offsetMinusOne = boost::integer_traits<T>::is_signed &&
                                boost::integer_traits<T>::const_min == i;

    if (offsetMinusOne)
      ++i;

    const bool negative = i<0;
    typename boost::make_unsigned<T>::type magnitude = std::abs(i);
    int loc = 0;
    while(magnitude>0)
    {
      packer.limbs()[loc] = magnitude & GMP_NUMB_MASK;
      magnitude /= GMP_NUMB_MAX;
      ++loc;
    }

    packer.commit();
    size = negate(negative, loc);

    if (offsetMinusOne)
      --(*this);
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

    const typename boost::make_unsigned<T>::type absI = std::abs(i);
    if ((!boost::integer_traits<T>::is_signed || boost::integer_traits<T>::const_min != i) &&
        absI <= GMP_NUMB_MAX)
    {
      // i is too small to be equal
      if (width() > 1)
        return false;

      // Both values fit into a single limb
      const mp_limb_t limb = static_cast<mp_limb_t>(absI);

      ConstPacker packer(this);
      return (packer.limbs()[0] & GMP_NUMB_MASK) == limb;
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

    const typename boost::make_unsigned<T>::type absI = std::abs(i);
    if ((!boost::integer_traits<T>::is_signed || boost::integer_traits<T>::const_min != i) &&
        absI <= GMP_NUMB_MAX)
    {
      // i has a smaller magnitude
      if (width() > 1)
        return isNegative();

      // Both values fit into a single limb
      ConstPacker packer(this);
      mp_limb_t left = packer.limbs()[0];
      mp_limb_t right = static_cast<mp_limb_t>(absI);

      if (isNegative())
        std::swap(left, right);

      return left < right;
    }

    return *this < Integer(i);
  }

  Integer(const int size, uintptr_t data);

public:
  void mpzInit(mpz_t& mpz) const;

  static bool isPrime(const Integer& x, int reps = 25);
  static Integer nextPrime(const Integer& x);
  static Integer gcd(const Integer& x, const Integer& y);
  static Integer lcm(const Integer& x, const Integer& y);

  Integer();
  ~Integer();
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
