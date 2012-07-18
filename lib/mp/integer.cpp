#include <gmp.h>
#include <apr_atomic.h>
#include <excafe/mp/integer.hpp>
#include <excafe/exception.hpp>
#include <excafe/numeric/cast.hpp>
#include <excafe/util/hash.hpp>
#include <excafe/util/hybrid_array.hpp>
#include <cstring>
#include <cassert>
#include <ostream>

namespace excafe
{

namespace mp
{

const Integer Integer::ONE =       Integer(1);
const Integer Integer::MINUS_ONE = Integer(-1);

void Integer::validate() const
{
  if (size != 0)
  {
    ConstPacker tp(this);
    assert(tp.limbs()[width()-1] != 0);
  }
}

Integer::Integer(const int _size, uintptr_t _data) :
  size(_size), data(_data)
{
  if ((data & VALUE_TAG) == 0)
    incrementUseCount(reinterpret_cast<mp_limb_t*>(data));

  EXCAFE_VALIDATE_INTEGER;
}

Integer::Integer() : size(0), data(VALUE_TAG)
{
  EXCAFE_VALIDATE_INTEGER;
}

Integer::Integer(const char* str) : size(0), data(VALUE_TAG)
{
  const std::size_t length = strlen(str);

  if (length == 0)
    CFD_EXCEPTION("Cannot construct integer from empty string.");

  const bool negative = (*str == '-');

  if (length == 1 && negative)
    CFD_EXCEPTION("Expected digits after minus sign when constructing integer.");

  for(std::size_t digit = (negative ? 1 : 0); digit < length; ++digit)
  {
    (*this) *= 10;
    const int digitValue = str[digit] - '0';

    if (digitValue < 0 || digitValue > 9)
      CFD_EXCEPTION("Invalid characters in integer string.");
    else
      (*this) += digitValue;
  }

  if (negative)
    (*this) = -(*this);

  EXCAFE_VALIDATE_INTEGER;
}

Integer::Integer(const Integer& i) : size(i.size), data(i.data)
{
  if ((data & VALUE_TAG) == 0)
    incrementUseCount(reinterpret_cast<mp_limb_t*>(data));

  EXCAFE_VALIDATE_INTEGER;
}

Integer::Integer(const int i) : size(0), data(VALUE_TAG)
{
  if (i != 0)
    initialise(i);

  EXCAFE_VALIDATE_INTEGER;
}

Integer::Integer(const long i) : size(0), data(VALUE_TAG)
{
  if (i != 0)
    initialise(i);

  EXCAFE_VALIDATE_INTEGER;
}

Integer::Integer(const unsigned int i) : size(0), data(VALUE_TAG)
{
  if (i != 0)
    initialise(i);

  EXCAFE_VALIDATE_INTEGER;
}

Integer::Integer(const unsigned long i) : size(0), data(VALUE_TAG)
{
  if (i != 0)
    initialise(i);

  EXCAFE_VALIDATE_INTEGER;
}

Integer::~Integer()
{
  if ((data & VALUE_TAG) == 0)
    decrementUseCount(reinterpret_cast<mp_limb_t*>(data));
}

void Integer::performAddition(const Integer* u, const Integer* v)
{
  int ausize = std::abs(u->size);
  int avsize = std::abs(v->size);

  // Ensure that u is always the larger operand
  if (ausize < avsize)
  {
    std::swap(u, v);
    std::swap(ausize, avsize);
  }

  Packer tp(this);
  const bool sameSign = (u->size < 0) == (v->size < 0);

  // up and vp must be constructed *after* the reallocation.
  if (sameSign)
  {
    // Possible carry
    const int maxWidth = ausize + 1;

    tp.reallocUnique(maxWidth);
    ConstPacker up(u);
    ConstPacker vp(v);

    mp_limb_t carry;
    if (avsize == 1)
      carry = mpn_add_1(tp.limbs(), up.limbs(), ausize, *vp.limbs());
    else
      carry = mpn_add(tp.limbs(), up.limbs(), ausize, vp.limbs(), avsize);

    tp.limbs()[ausize] = carry;
    size = negate(u->size < 0, ausize + carry);
  }
  else
  {
    tp.reallocUnique(ausize);
    ConstPacker up(u);
    ConstPacker vp(v);

    if (ausize != avsize)
    {
      mpn_sub(tp.limbs(), up.limbs(), ausize, vp.limbs(), avsize);
      size = negate(u->size < 0, tp.computeWidth(ausize));
    }
    else if (mpn_cmp(up.limbs(), vp.limbs(), ausize) < 0)
    {
      mpn_sub_n(tp.limbs(), vp.limbs(), up.limbs(), ausize);
      size = negate(v->size < 0, tp.computeWidth(ausize));
    }
    else
    {
      mpn_sub_n(tp.limbs(), up.limbs(), vp.limbs(), ausize);
      size = negate(u->size < 0, tp.computeWidth(ausize));
    }
  }

  tp.commit();
}

std::size_t Integer::numLimbs(const int bits)
{
  return (bits + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
}

mp_limb_t* Integer::allocateLimbs(const std::size_t count)
{
  void* data = malloc(sizeof(Header) + sizeof(mp_limb_t)*count);
  static_cast<Header*>(data)->allocated = count;
  apr_atomic_set32(&static_cast<Header*>(data)->count, 1);
  return static_cast<mp_limb_t*>(static_cast<void*>(static_cast<char*>(data)+sizeof(Header)));
}

std::size_t Integer::countLimbs(const mp_limb_t* limbs)
{
  const void* data = static_cast<const char*>(static_cast<const void*>(limbs)) - sizeof(Header);
  return static_cast<const Header*>(data)->allocated;
}

void Integer::incrementUseCount(mp_limb_t* limbs)
{
  void* data = static_cast<char*>(static_cast<void*>(limbs)) - sizeof(Header);
  apr_atomic_inc32(&static_cast<Header*>(data)->count);
}

void Integer::decrementUseCount(mp_limb_t* limbs)
{
  void* data = static_cast<char*>(static_cast<void*>(limbs)) - sizeof(Header);
  const int result = apr_atomic_dec32(&static_cast<Header*>(data)->count);

  if (result == 0)
    free(data);
}

bool Integer::isUnique(mp_limb_t* limbs)
{
  void* data = static_cast<char*>(static_cast<void*>(limbs)) - sizeof(Header);
  const int result = apr_atomic_read32(&static_cast<Header*>(data)->count);

  return result == 1;
}

void Integer::mpzInit(mpz_t& mpz) const
{
  ConstPacker tp(this);
  const int w = width();

  // Add one for the sign bit
  mpz_init2(mpz, w * GMP_NUMB_BITS + 1);

  for(int i=0; i<w; ++i)
  {
    mpz_mul_2exp(mpz, mpz, GMP_NUMB_BITS);
    mpz_add_ui(mpz, mpz, tp.limbs()[i]);
  }
  
  if (size < 0)
    mpz_neg(mpz, mpz);
}

Integer Integer::gcd(const Integer& x, const Integer& y)
{
  if (x.size == 0)
    return y;
  else if (y.size == 0)
    return x;

  Integer result;
  Packer rp(&result);
  ConstPacker xp(&x);
  ConstPacker yp(&y);

  if (y.width() == 1)
  {
    rp.reallocUnique(1);
    result.size = 1;
    rp.limbs()[0] = mpn_gcd_1(xp.limbs(), x.width(), yp.limbs()[0]);
  }
  else if (x.width() == 1)
  {
    rp.reallocUnique(1);
    result.size = 1;
    rp.limbs()[0] = mpn_gcd_1(yp.limbs(), y.width(), xp.limbs()[0]);
  }
  else
  {
    Integer cx(x), cy(y);

    // Make x the largest
    if (cx.width() < cy.width())
      cx.swap(cy);

    // mpn_gcd clobbers
    Packer cxp(&cx);
    Packer cyp(&cy);
    cxp.reallocUnique(cx.width());
    cyp.reallocUnique(cy.width());

    const mp_bitcnt_t cxZeros = mpn_scan1(cxp.limbs(), 0);
    const mp_bitcnt_t cyZeros = mpn_scan1(cyp.limbs(), 0);

    // Remove leading zeros from cy
    cy >>= cyZeros;
    rp.reallocUnique(cy.width());
    result.size = cy.width();
    cyp.unpack();
    mpn_gcd(rp.limbs(), cxp.limbs(), cx.width(), cyp.limbs(), cy.width());
    cyp.abort();
    cxp.abort();

    // Multiply result by common powers of two between cx and cy
    result <<= std::min(cxZeros, cyZeros);
  }

  rp.commit();
  return result;
}

Integer Integer::lcm(const Integer& x, const Integer& y)
{
  return x*(y/gcd(x, y));
}

bool Integer::operator==(const Integer& i) const
{
  if (size == i.size)
  {
    if (size == 0)
      return true;

    ConstPacker tp(this);
    ConstPacker ip(&i);
    return mpn_cmp(tp.limbs(), ip.limbs(), width()) == 0;
  }
  return false;
}

bool Integer::operator==(const int i) const
{
  return equal(i);
}

bool Integer::operator==(const long i) const
{
  return equal(i);
}

bool Integer::operator<(const Integer& i) const
{
  if (size == i.size)
  {
    if (size == 0)
    {
      return false;
    }
    else if (isNegative() != i.isNegative())
    {
      return isNegative();
    }
    else
    {
      const int comparison = (isNegative() ? 1 : -1);
      ConstPacker tp(this);
      ConstPacker ip(&i);
      return mpn_cmp(tp.limbs(), ip.limbs(), width()) == comparison;
    }
  }
  else
  {
    return size < i.size;
  }
}

bool Integer::operator<(const int i) const
{
  return lessThan(i);
}

bool Integer::operator<(const long i) const
{
  return lessThan(i);
}

Integer& Integer::operator=(const Integer& i)
{
  if ((i.data & VALUE_TAG) == 0)
    incrementUseCount(reinterpret_cast<mp_limb_t*>(i.data));

  if ((data & VALUE_TAG) == 0)
    decrementUseCount(reinterpret_cast<mp_limb_t*>(data));

  size = i.size;
  data = i.data;

  EXCAFE_VALIDATE_INTEGER;
  return *this;
}

Integer& Integer::operator+=(const Integer& i)
{
  if (i.size == 0)
    return *this;
  else if (size == 0)
    return *this = i;
  else
    performAddition(this, &i);

  EXCAFE_VALIDATE_INTEGER;
  return *this;
}

Integer& Integer::operator-=(const Integer& i)
{
  if (i.size == 0)
    return *this;
  else if (size == 0)
    return *this = -i;
  else
  {
    const Integer negated = -i;
    performAddition(this, &negated);
  }

  EXCAFE_VALIDATE_INTEGER;
  return *this;
}

Integer& Integer::operator++()
{
  if (size == 0)
    *this = ONE;
  else
    performAddition(this, &ONE);

  EXCAFE_VALIDATE_INTEGER;
  return *this;
}

Integer& Integer::operator--()
{
  if (size == 0)
    *this = MINUS_ONE;
  else
    performAddition(this, &MINUS_ONE);

  EXCAFE_VALIDATE_INTEGER;
  return *this;
}

Integer& Integer::operator*=(const Integer& i)
{
  if (size == 0 || i.size == 0)
  {
    size = 0;
  }
  else
  {
    const bool negative = isNegative() ^ i.isNegative();

    const Integer* u = this;
    const Integer* v = &i;

    int ausize = u->width();
    int avsize = v->width();
    const int maxWidth = ausize + avsize;

    // Ensure that u is always the larger operand
    if (ausize < avsize)
    {
      std::swap(u, v);
      std::swap(ausize, avsize);
    }

    // mpn_mul does not allow any operand to overlap with the result
    Integer result;
    Packer rp(&result);
    ConstPacker up(u);
    ConstPacker vp(v);

    rp.reallocUnique(maxWidth);

    if (avsize == 1)
      rp.limbs()[maxWidth-1] = mpn_mul_1(rp.limbs(), up.limbs(), ausize, *vp.limbs());
    else
      mpn_mul(rp.limbs(), up.limbs(), ausize, vp.limbs(), avsize);

    result.size = negate(negative, rp.computeWidth(maxWidth));
    rp.commit();
  
    swap(result);
  }

  EXCAFE_VALIDATE_INTEGER;
  return *this;
}

Integer& Integer::operator/=(const Integer& dividend)
{
  if (size == 0 || width() < dividend.width())
  {
    size = 0;
  }
  else
  {
    const bool negative = isNegative() ^ dividend.isNegative();
    const int maxWidth = width() - dividend.width() + 1;
    ConstPacker dp(&dividend);
    Packer tp(this);

    // Ensure that we are unique since we intend to overwite our data
    tp.reallocUnique(width());

    // Perform the division, overwriting our original value with the remainder
    if (dividend.width() == 1)
    {
      mpn_divrem_1(tp.limbs(), 0, tp.limbs(), width(), *dp.limbs());
    }
    else
    {
      // Allocate space for our result
      util::HybridArray<mp_limb_t, STACK_LIMBS> result(maxWidth);

      mpn_tdiv_qr(result.get(), tp.limbs(), 0, tp.limbs(), width(), dp.limbs(), dividend.width());

      // Copy result of division back to data
      memcpy(tp.limbs(), result.get(), maxWidth*sizeof(mp_limb_t));
    }

    size = negate(negative, tp.computeWidth(maxWidth));
    tp.commit();
  }

  EXCAFE_VALIDATE_INTEGER;
  return *this;
}

Integer& Integer::operator%=(const Integer& dividend)
{
  if (size != 0 || width() >= dividend.width())
  {
    const bool negative = isNegative();
    const int maxWidth = width() - dividend.width() + 1;
    ConstPacker dp(&dividend);
    Packer tp(this);

    // Ensure that we are unique since we intend to overwite our data
    tp.reallocUnique(width());

    if (dividend.width() == 1)
    {
      tp.limbs()[0] = mpn_mod_1(tp.limbs(), width(), *dp.limbs());
    }
    else
    {
      // Allocate space for our result
      util::HybridArray<mp_limb_t, STACK_LIMBS> result(maxWidth);

      // Perform the division, overwriting our original value with the remainder
      mpn_tdiv_qr(result.get(), tp.limbs(), 0, tp.limbs(), width(), dp.limbs(), dividend.width());
    }

    // Compute size of remainder
    size = negate(negative, tp.computeWidth(dividend.width()));

    // Commit changes
    tp.commit();
  }

  EXCAFE_VALIDATE_INTEGER;
  return *this;
}

Integer& Integer::operator>>=(unsigned i)
{
  const int oldWidth = width();
  const unsigned limbCount = i/GMP_NUMB_BITS;

  if (limbCount >= oldWidth)
  {
    size = 0;
    return *this;
  }

  i -= GMP_NUMB_BITS*limbCount;
  const int newWidth = oldWidth - limbCount;
  Packer tp(this);

  // Make unique
  tp.reallocUnique(oldWidth);

  if (i == 0)
    mpn_copyi(tp.limbs(), tp.limbs()+limbCount, newWidth);
  else
    mpn_rshift(tp.limbs(), tp.limbs()+limbCount, newWidth, i);

  size = negate(isNegative(), tp.computeWidth(newWidth));
  tp.commit();

  EXCAFE_VALIDATE_INTEGER;
  return *this;
}

Integer& Integer::operator<<=(unsigned i)
{
  const int oldWidth = width();
  const unsigned limbCount = i/GMP_NUMB_BITS;
  i -= GMP_NUMB_BITS*limbCount;
  const int newWidth = oldWidth + limbCount + (i>0 ? 1 : 0);
  Packer tp(this);

  // Enlarge and make unique
  tp.reallocUnique(newWidth);

  if (i == 0)
    mpn_copyd(tp.limbs()+limbCount, tp.limbs(), oldWidth);
  else
    tp.limbs()[newWidth-1] = mpn_lshift(tp.limbs()+limbCount, tp.limbs(), oldWidth, i);

  mpn_zero(tp.limbs(), limbCount);

  size = negate(isNegative(), tp.computeWidth(newWidth));
  tp.commit();

  EXCAFE_VALIDATE_INTEGER;
  return *this;
}

Integer Integer::operator-() const
{
  return Integer(-size, data);
}

Integer Integer::abs() const
{
  return Integer(width(), data);
}

Integer Integer::isqrt() const
{
  if (size < 0)
    CFD_EXCEPTION("Can only take sqrt of positive integer.");

  const int resultWidth = width()/2 + width()%2;
  Integer result;
  Packer rp(&result);
  ConstPacker tp(this);

  rp.reallocUnique(resultWidth);
  mpn_sqrtrem(rp.limbs(), NULL, tp.limbs(), width());
  result.size = rp.computeWidth(resultWidth);
  rp.commit();

  return result;
}

Integer Integer::pow(const long exponent) const
{
  // Exponentiation by squaring
  switch(exponent)
  {
    case 0: return Integer(1);
    case 1: return *this;
    case 2: return (*this)*(*this);
    default:
    {
      if (exponent < 0)
      {
        CFD_EXCEPTION("Cannot raise integer to negative exponent.");
      }
      else if (exponent%2 == 1)
      {
        const Integer other = pow((exponent-1)/2);
        return (*this) * other * other;
      }
      else
      {
        const Integer subTerm = pow(exponent/2);
        return subTerm * subTerm;
      }
    }
  }
}

std::size_t Integer::hash() const
{
  std::size_t result = 0x42cbce0c;
  util::hash_accum(result, size < 0);

  ConstPacker tp(this);
  const std::size_t w = width();
  for(std::size_t i=0; i<w; ++i)
    util::hash_accum(result, tp.limbs()[i] & GMP_NUMB_MASK);

  return result;
}

void Integer::write(std::ostream& out) const
{
  if (size == 0)
  {
    out << "0";
  }
  else
  {
    const int base = 10;
    const int bitsPerCharFloor = 3;
    const long bits = width()*sizeof(mp_limb_t)*CHAR_BIT;
    const long maxChars = bits/bitsPerCharFloor + 1;

    if (size < 0)
      out << "-";

    Integer clobbered(*this);
    Packer cp(&clobbered);
    cp.reallocUnique(width());

    std::vector<unsigned char> characters(maxChars);
    const mp_size_t chars = mpn_get_str(&characters[0], base, cp.limbs(), clobbered.width());
    cp.abort();

    mp_size_t loc = 0;
    // Skip any leading zeros
    while (loc < chars && characters[loc] == 0) ++loc;

    while(loc < chars)
      out.put('0' + characters[loc++]);
  }
}

void Integer::swap(Integer& i)
{
  std::swap(size, i.size);
  std::swap(data, i.data);
}

long Integer::toLong() const
{
  mpz_t mpz;
  mpzInit(mpz);
  const long value = mpz_get_ui(mpz);
  mpz_clear(mpz);

  if (mpz_fits_slong_p(mpz))
    return size < 0 ? -value : value;
  else
    CFD_EXCEPTION("Cannot represent Integer as long.");
}

int Integer::toInt() const
{
  return excafe::numeric_cast<int>(toLong());
}

double Integer::toDouble() const
{
  mpz_t mpz;
  mpzInit(mpz);
  const double result = mpz_get_d(mpz);
  mpz_clear(mpz);
  return result;
}

float Integer::toFloat() const
{
  return excafe::numeric_cast<float>(toDouble());
}

std::size_t hash_value(const Integer& i)
{
  return i.hash();
}

Integer abs(const Integer& i)
{
  return i.abs();
}

Integer pow(const Integer& i, const int exponent)
{
  return i.pow(exponent);
}

Integer pow(const Integer& i, const long exponent)
{
  return i.pow(exponent);
}

Integer isqrt(const Integer& i)
{
  return i.isqrt();
}

std::ostream& operator<<(std::ostream& o, const Integer& i)
{
  i.write(o);
  return o;
}

}

}
