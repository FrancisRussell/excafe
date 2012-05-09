// We include this before gmp.h to enable gmp's printf functionality
#include <cstdarg>

#include <gmp.h>
#include <cstring>
#include <ostream>
#include <excafe/exception.hpp>
#include <excafe/mp/float.hpp>
#include <excafe/mp/integer.hpp>
#include <excafe/mp/rational.hpp>
#include <excafe/numeric/cast.hpp>
#include <excafe/util/hash.hpp>

namespace excafe
{

namespace mp
{

Float::Float()
{
  mpf_init(value);
}

Float::Float(const char* str)
{
  const int result = mpf_init_set_str(value, str, 10);

  if (result != 0)
    CFD_EXCEPTION("Invalid specification for floating point value.");
}

Float::Float(const double f)
{
  mpf_init_set_d(value, f);
}

Float::Float(const Float& f)
{
  mpf_init_set(value, f.value);
}

Float::Float(const Integer& i)
{
  mpf_init(value);

  mpz_t mpz;
  i.mpzInit(mpz);
  mpf_set_z(value, mpz);
  mpz_clear(mpz);
}

Float::Float(const Rational& r)
{
  mpf_init(value);

  mpq_t mpq;
  r.mpqInit(mpq);
  mpf_set_q(value, mpq);
  mpq_clear(mpq);
}

Float& Float::operator=(const Float& f)
{
  mpf_set(value, f.value);
  return *this;
}

bool Float::operator==(const Float& f) const
{
  return mpf_cmp(value, f.value) == 0;
}

bool Float::operator<(const Float& f) const
{  
  return mpf_cmp(value, f.value) < 0;
}

Float& Float::operator+=(const Float& f)
{
  mpf_add(value, value, f.value);
  return *this;
}

Float& Float::operator-=(const Float& f)
{
  mpf_sub(value, value, f.value);
  return *this;
}

Float& Float::operator*=(const Float& f)
{
  mpf_mul(value, value, f.value);
  return *this;
}

Float& Float::operator/=(const Float& f)
{
  mpf_div(value, value, f.value);
  return *this;
}

Float Float::operator-() const
{
  Float result;
  mpf_neg(result.value, value);
  return result;
}

Float Float::abs() const
{
  Float result;
  mpf_abs(result.value, value);
  return result;
}

Float Float::pow(const long exponent) const
{
  Float result(*this);

  if (exponent < 0)
    mpf_ui_div(result.value, 1, result.value);

  mpf_pow_ui(result.value, result.value, std::abs(exponent));
  return result;
}

void Float::write(std::ostream& out) const
{
  char* buffer;
  const int bufferSize = gmp_asprintf(&buffer, "%.Fe", value);

  out << buffer;

  void (*freefunc)(void*, size_t);
  mp_get_memory_functions (NULL, NULL, &freefunc);
  freefunc(buffer, bufferSize);
}

float Float::toFloat() const
{
  return excafe::numeric_cast<float>(toDouble());
}

double Float::toDouble() const
{
  return mpf_get_d(value);
}

Float::~Float()
{
  mpf_clear(value);
}

std::ostream& operator<<(std::ostream& o, const Float& f)
{
  f.write(o);
  return o;
}

Float abs(const Float& f)
{
  return f.abs();
}

Float pow(const Float& f, const long exponent)
{
  return f.pow(exponent);
}

Float pow(const Float& f, const int exponent)
{
  return f.pow(exponent);
}

void Float::swap(Float& f)
{
  mpf_swap(value, f.value);
}

}

}
