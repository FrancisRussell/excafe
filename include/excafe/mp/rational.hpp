#ifndef SIMPLE_CFD_MP_RATIONAL_HPP
#define SIMPLE_CFD_MP_RATIONAL_HPP

#include "mp_fwd.hpp"
#include "integer.hpp"
#include <boost/operators.hpp>

#ifdef EXCAFE_VALIDATE_MP
#define EXCAFE_VALIDATE_RATIONAL this->validate();
#else
#define EXCAFE_VALIDATE_RATIONAL
#endif


namespace excafe
{

namespace mp
{

class Rational : boost::totally_ordered<Rational>,
                 boost::arithmetic<Rational>
{
private:
  friend class Integer;
  friend class Float;

  Integer numerator;
  Integer denominator;

  struct no_normalise_tag {};

  Rational(const Integer& _numerator, const Integer& _denominator, struct no_normalise_tag);
  void normalise();
  void validate() const;
  void mpqInit(mpq_t& mpq) const;

public:
  Rational();
  Rational(const int i);
  Rational(const long i);
  Rational(const unsigned int i);
  Rational(const unsigned long i);
  Rational(const Integer& i);
  Rational(const int numerator, const int denominator);
  Rational(const long numerator, const long denominator); 
  Rational(const Integer& numerator, const Integer& denominator);

  Integer getNumerator() const;
  Integer getDenominator() const;

  bool operator==(const Rational& r) const;
  bool operator==(const Integer& i) const;
  bool operator==(int i) const;
  bool operator==(long i) const;
  bool operator<(const Rational& r) const;
  bool operator<(const Integer& r) const;
  bool operator<(int i) const;
  bool operator<(long i) const;
  Rational operator-() const;
  Rational& operator+=(const Rational& r);
  Rational& operator-=(const Rational& r);
  Rational& operator*=(const Rational& r);
  Rational& operator/=(const Rational& r);
  Rational& operator++();
  Rational& operator--();

  Rational reciprocal() const;
  Rational abs() const;
  Rational pow(const long exponent) const;
  void write(std::ostream& out) const;
  std::size_t hash() const;
  double toDouble() const;
  float toFloat() const;
};

std::ostream& operator<<(std::ostream& o, const Rational& r);
std::size_t hash_value(const Rational& r);
Rational abs(const Rational& r);
Rational pow(const Rational& r, const int exponent);
Rational pow(const Rational& r, const long exponent);

}

}

#endif
