#include <gmp.h>
#include <simple_cfd/mp/integer.hpp>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/numeric/cast.hpp>
#include <simple_cfd/util/hash.hpp>
#include <cstring>
#include <cassert>
#include <ostream>

namespace cfd
{

namespace mp
{

Integer::Integer(const int _size, const int _allocated, const boost::shared_array<mp_limb_t>& _data) :
  size(_size), allocated(_allocated), data(_data)
{
}

Integer::Integer() : size(0), allocated(0)
{
}

Integer::Integer(const char* str) : size(0), allocated(0)
{
  const std::size_t length = strlen(str);

  if (length == 0)
    return;

  const bool negative = (*str == '-');

  for(std::size_t digit = (negative ? 1 : 0); digit < length; ++digit)
  {
    (*this) *= 10;
    const int digitValue = str[digit] - '0';

    if (digitValue < 0 || digitValue > 9)
      CFD_EXCEPTION("Invalid characters in integer string");
    else
      (*this) += digitValue;
  }

  if (negative)
    (*this) = -(*this);
}

Integer::Integer(const Integer& i) : size(i.size), allocated(i.allocated), data(i.data)
{
}

Integer::Integer(const int i) : size(0), allocated(0)
{
  if (i != 0)
    initialise(i);
}

Integer::Integer(const long i) : size(0), allocated(0)
{
  if (i != 0)
    initialise(i);
}

Integer::Integer(const unsigned int i) : size(0), allocated(0)
{
  if (i != 0)
    initialise(i);
}

Integer::Integer(const unsigned long i) : size(0), allocated(0)
{
  if (i != 0)
    initialise(i);
}

void Integer::reallocUnique(const int newAllocated)
{
  if (newAllocated > allocated || !data.unique())
  {
    mp_limb_t* const newLimbs = new mp_limb_t[newAllocated];
    memcpy(newLimbs, limbs(), std::min(allocated, newAllocated) * sizeof(mp_limb_t));
    allocated = newAllocated;
    data.reset(newLimbs);
  }
}

int Integer::computeWidth(const int maxWidth) const
{
  int cwidth;
  for(cwidth = maxWidth; cwidth>0; --cwidth)
  {
    if (data[cwidth-1] != 0)
      break;
  }
  return cwidth;
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

  const bool sameSign = (u->size < 0) == (v->size < 0);
  if (sameSign)
  {
    // Possible carry
    const int maxWidth = ausize + 1;
    reallocUnique(maxWidth);

    const mp_limb_t carry = mpn_add(limbs(), u->limbs(), ausize, v->limbs(), avsize);
    limbs()[ausize] = carry;
    size = negate(u->size < 0, ausize + carry);
  }
  else
  {
    reallocUnique(ausize);

    if (ausize != avsize)
    {
      mpn_sub(limbs(), u->limbs(), ausize, v->limbs(), avsize);
      size = negate(u->size < 0, computeWidth(ausize));
    }
    else if (mpn_cmp(u->limbs(), v->limbs(), ausize) < 0)
    {
      mpn_sub_n(limbs(), v->limbs(), u->limbs(), ausize);
      size = negate(v->size < 0, computeWidth(ausize));
    }
    else
    {
      mpn_sub_n(limbs(), u->limbs(), v->limbs(), ausize);
      size = negate(u->size < 0, computeWidth(ausize));
    }
  }
}

std::size_t Integer::numLimbs(const int bits)
{
  std::size_t result = bits / GMP_NUMB_BITS;
  result += ((bits % GMP_NUMB_BITS) != 0);
  return result;
}

void Integer::mpzInit(mpz_t& mpz) const
{
  const int w = width();

  // Add one for the sign bit
  mpz_init2(mpz, w * GMP_NUMB_BITS + 1);

  for(int i=0; i<w; ++i)
  {
    mpz_mul_2exp(mpz, mpz, GMP_NUMB_BITS);
    mpz_add_ui(mpz, mpz, limbs()[i]);
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
  if (y.width() == 1)
  {
    result.reallocUnique(1);
    result.size = 1;
    *result.limbs() = mpn_gcd_1(x.limbs(), x.width(), *y.limbs());
  }
  else if (x.width() == 1)
  {
    result.reallocUnique(1);
    result.size = 1;
    *result.limbs() = mpn_gcd_1(y.limbs(), y.width(), *x.limbs());
  }
  else
  {
    Integer cx(x), cy(y);

    // Make x the largest
    if (cx.width() < cy.width())
      cx.swap(cy);

    // mpn_gcd clobbers
    cx.reallocUnique(cx.width());
    cy.reallocUnique(cy.width());

    const mp_bitcnt_t cxZeros = mpn_scan1(cx.limbs(), 0);
    const mp_bitcnt_t cyZeros = mpn_scan1(cy.limbs(), 0);

    // Remove leading zeros from cy
    cy >>= cyZeros;
    result.reallocUnique(cy.width());
    result.size = cy.width();
    mpn_gcd(result.limbs(), cx.limbs(), cx.width(), cy.limbs(), cy.width());

    // Multiply result by common powers of two between cx and cy
    result <<= std::min(cxZeros, cyZeros);
  }

  assert(result != 0);
  return result;
}

Integer Integer::lcm(const Integer& x, const Integer& y)
{
  return x*(y/gcd(x, y));
}

bool Integer::operator==(const Integer& i) const
{
  return size == i.size 
         && (size == 0 || mpn_cmp(limbs(), i.limbs(), width()) == 0);
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
    return size != 0
           && (isNegative() ^ (mpn_cmp(limbs(), i.limbs(), width()) < 0));
  else
    return size < i.size;
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
  size = i.size;
  allocated = i.allocated;
  data = i.data;

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
    const Integer negated = -(*this);
    performAddition(this, &negated);
  }

  return *this;
}

Integer& Integer::operator++()
{
  if (size == 0)
  {
    initialise(1);
  }
  else
  {
    const Integer one(1);
    performAddition(this, &one);
  }

  return *this;
}

Integer& Integer::operator--()
{
  if (size == 0)
  {
    initialise(-1);
  }
  else
  {
    const Integer minusOne(-1);
    performAddition(this, &minusOne);
  }

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
    const int maxWidth = width() + i.width();
    reallocUnique(maxWidth);
    mpn_mul(limbs(), limbs(), width(), i.limbs(), i.width());
    size = negate(negative, computeWidth(maxWidth));
  }

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

    // Ensure that we are unique since we intend to overwite our data
    reallocUnique(width());

    // Allocate space for our result
    util::HybridArray<mp_limb_t, STACK_LIMBS> result(maxWidth);

    // Perform the division, overwriting our original value with the remainder
    mpn_tdiv_qr(result.get(), limbs(), 0, limbs(), width(), dividend.limbs(), dividend.width());

    // Copy result of division back to data
    memcpy(limbs(), result.get(), maxWidth*sizeof(mp_limb_t));
    size = negate(negative, computeWidth(maxWidth));
  }

  return *this;
}

Integer& Integer::operator%=(const Integer& dividend)
{
  if (size != 0 || width() >= dividend.width())
  {
    const bool negative = isNegative();
    const int maxWidth = width() - dividend.width() + 1;

    // Ensure that we are unique since we intend to overwite our data
    reallocUnique(width());

    // Allocate space for our result
    util::HybridArray<mp_limb_t, STACK_LIMBS> result(maxWidth);

    // Perform the division, overwriting our original value with the remainder
    mpn_tdiv_qr(result.get(), limbs(), 0, limbs(), width(), dividend.limbs(), dividend.width());

    // Compute size of remainder
    size = negate(negative, computeWidth(dividend.width()));
  }

  return *this;
}

Integer& Integer::operator>>=(unsigned i)
{
  const int oldWidth = width();

  // Make unique
  reallocUnique(oldWidth);

  while(i>0)
  {
    const unsigned shift = std::min(i, static_cast<unsigned>(GMP_NUMB_BITS));
    mpn_rshift(limbs(), limbs(), width(), shift);
    i -= shift;
  }

  size = negate(isNegative(), computeWidth(oldWidth - i/GMP_NUMB_BITS));
  return *this;
}

Integer& Integer::operator<<=(unsigned i)
{
  const int newWidth = width() + i/GMP_NUMB_BITS + 1;

  // Enlarge and make unique
  reallocUnique(newWidth);

  while(i>0)
  {
    const unsigned shift = std::min(i, static_cast<unsigned>(GMP_NUMB_BITS));
    limbs()[width()] = mpn_lshift(limbs(), limbs(), width(), shift);
    size += (size < 0 ? -1 : 1);
    i -= shift;
  }

  size = negate(isNegative(), computeWidth(newWidth));
  return *this;
}

Integer Integer::operator-() const
{
  return Integer(-size, allocated, data);
}

Integer Integer::abs() const
{
  return Integer(width(), allocated, data);
}

Integer Integer::isqrt() const
{
  if (size < 0)
    CFD_EXCEPTION("Can only take sqrt of positive integer.");

  Integer result;
  result.reallocUnique(width()/2 + 1);
  result.size = mpn_sqrtrem(result.limbs(), NULL, limbs(), width());
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

  const std::size_t w = width();
  for(std::size_t i=0; i<w; ++i)
    util::hash_accum(result, data[i] & GMP_NUMB_MASK);

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
    clobbered.reallocUnique(width());

    std::vector<unsigned char> characters(maxChars);
    const mp_size_t chars = mpn_get_str(&characters[0], base, clobbered.limbs(), clobbered.width());

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
  std::swap(allocated, i.allocated);
  data.swap(i.data);
}

long Integer::toLong() const
{
  mpz_t mpz;
  mpzInit(mpz);
  const long value = mpz_get_ui(mpz);
  mpz_clear(mpz);

  if (value > 10000)
  {
    std::cout << "hmm: " << *this << std::endl;
    assert(false);
  }

  if (mpz_fits_slong_p(mpz))
    return size < 0 ? -value : value;
  else
    CFD_EXCEPTION("Cannot represent Integer as long.");
}

int Integer::toInt() const
{
  return cfd::numeric_cast<int>(toLong());
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
  return cfd::numeric_cast<float>(toDouble());
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
