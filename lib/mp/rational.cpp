#include <excafe/mp/rational.hpp>
#include <excafe/mp/integer.hpp>
#include <excafe/numeric/cast.hpp>
#include <excafe/util/hash.hpp>
#include <gmp.h>
#include <ostream>
#include <cassert>

namespace excafe
{

namespace mp
{
  
Rational::Rational(const Integer& _numerator, const Integer& _denominator, struct no_normalise_tag) :
  numerator(_numerator), denominator(_denominator)
{
  EXCAFE_VALIDATE_RATIONAL
}

Rational::Rational() : numerator(0), denominator(1)
{
  EXCAFE_VALIDATE_RATIONAL
}

Rational::Rational(const int i) : numerator(i), denominator(1)
{
  EXCAFE_VALIDATE_RATIONAL
}

Rational::Rational(const long i) : numerator(i), denominator(1)
{
  EXCAFE_VALIDATE_RATIONAL
}

Rational::Rational(const unsigned int i) : numerator(i), denominator(1)
{
  EXCAFE_VALIDATE_RATIONAL
}

Rational::Rational(const unsigned long i) : numerator(i), denominator(1)
{
  EXCAFE_VALIDATE_RATIONAL
}

Rational::Rational(const Integer& i) : numerator(i), denominator(1)
{
  EXCAFE_VALIDATE_RATIONAL
}

Rational::Rational(const int _numerator, const int _denominator) : 
  numerator(_numerator), denominator(_denominator)
{
  normalise();
  EXCAFE_VALIDATE_RATIONAL
}

Rational::Rational(const long _numerator, const long _denominator) : 
  numerator(_numerator), denominator(_denominator)
{
  normalise();
  EXCAFE_VALIDATE_RATIONAL
}

Rational::Rational(const Integer& _numerator, const Integer& _denominator) : 
  numerator(_numerator), denominator(_denominator)
{
  normalise();
  EXCAFE_VALIDATE_RATIONAL
}

void Rational::normalise()
{
  assert(denominator != 0);

  if (numerator == 0)
    denominator = 1;
  
  if (denominator == 1)
    return;

  if (denominator < 0)
  {
    numerator = -numerator;
    denominator = -denominator;
  }

  const Integer factor = Integer::gcd(numerator, denominator);
  if (factor != 1)
  {
    numerator /= factor;
    denominator /= factor;
  }
}

void Rational::validate() const
{
  Rational normalised(*this);
  normalised.normalise();

  assert(normalised == *this);
}

void Rational::mpqInit(mpq_t& mpq) const
{
  mpz_t num, den;
  numerator.mpzInit(num);
  denominator.mpzInit(den);

  mpq_init(mpq);
  mpq_set_num(mpq, num);
  mpq_set_den(mpq, den);

  mpz_clear(num);
  mpz_clear(den);
}

bool Rational::operator==(const Rational& r) const
{
  return numerator == r.numerator && denominator == r.denominator;
}

bool Rational::operator<(const Rational& r) const
{
  return numerator*r.denominator < r.numerator*denominator;
}

Integer Rational::getNumerator() const
{
  return numerator;
}

Integer Rational::getDenominator() const
{
  return denominator;
}

Rational Rational::operator-() const
{
  return Rational(-numerator, denominator, no_normalise_tag());
}

Rational& Rational::operator+=(const Rational& r)
{
  const Integer lcd = Integer::lcm(denominator, r.denominator);
  numerator = numerator*(lcd/denominator) + r.numerator*(lcd/r.denominator);
  denominator = lcd;

  normalise();

  EXCAFE_VALIDATE_RATIONAL
  return *this;
}

Rational& Rational::operator-=(const Rational& r)
{
  (*this) += -r;

  EXCAFE_VALIDATE_RATIONAL
  return *this;
}

Rational& Rational::operator*=(const Rational& r)
{
  const Integer gcd1 = Integer::gcd(numerator, r.denominator);
  const Integer gcd2 = Integer::gcd(r.numerator, denominator);

  numerator = (numerator/gcd1)*(r.numerator/gcd2);
  denominator = (denominator/gcd2)*(r.denominator/gcd1);

  EXCAFE_VALIDATE_RATIONAL
  return *this;
}

Rational Rational::reciprocal() const
{
  assert(numerator != 0);

  if (numerator < 0)
    return Rational(-denominator, -numerator, no_normalise_tag());
  else
    return Rational(denominator, numerator, no_normalise_tag());
}

Rational& Rational::operator/=(const Rational& r)
{
  (*this) *= r.reciprocal();

  EXCAFE_VALIDATE_RATIONAL
  return *this;
}

Rational& Rational::operator++()
{
  numerator += denominator;

  EXCAFE_VALIDATE_RATIONAL
  return *this;
}

Rational& Rational::operator--()
{
  numerator -= denominator;

  EXCAFE_VALIDATE_RATIONAL
  return *this;
}

Rational Rational::abs() const
{
  if (numerator < 0)
    return Rational(-numerator, denominator, no_normalise_tag());
  else
    return *this;
}

Rational Rational::pow(const long exponent) const
{
  const long absExp = std::abs(exponent);

  const Rational value = (exponent < 0 ? reciprocal() : *this);
  const Integer num = value.numerator.pow(absExp);
  const Integer den = value.denominator.pow(absExp);

  return Rational(num, den, no_normalise_tag());
}

void Rational::write(std::ostream& out) const
{
  out << numerator << "/" << denominator;
}

std::size_t Rational::hash() const
{
  std::size_t result = 0xde7fa5b9;
  util::hash_accum(result, numerator);
  util::hash_accum(result, denominator);
  return result;
}

double Rational::toDouble() const
{
  mpq_t mpq;
  mpqInit(mpq);
  const double result = mpq_get_d(mpq);
  mpq_clear(mpq);
  return result;
}

float Rational::toFloat() const
{
  return excafe::numeric_cast<float>(toDouble());
}

std::size_t hash_value(const Rational& r)
{
  return r.hash();
}

Rational abs(const Rational& r)
{
  return r.abs();
}

Rational pow(const Rational& r, const int exponent)
{
  return r.pow(exponent);
}

Rational pow(const Rational& r, const long exponent)
{
  return r.pow(exponent);
}

std::ostream& operator<<(std::ostream& o, const Rational& r)
{
  r.write(o);
  return o;
}

}

}
