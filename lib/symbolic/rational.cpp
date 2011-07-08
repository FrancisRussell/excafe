#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/make_expr_from.hpp>
#include <simple_cfd/numeric/math_utilities.hpp>
#include <simple_cfd/util/hash.hpp>
#include <ostream>
#include <cassert>
#include <cmath>
#include <climits>
#include <set>
#include <cln/integer.h>
#include <cln/rational.h>

namespace cfd
{

namespace symbolic
{

const Expr Rational::zero()
{
  static const Expr value = make_expr_from(Rational(0));
  return value;
}

const Expr Rational::one()
{
  static const Expr value = make_expr_from(Rational(1));
  return value;
}

Rational::Rational() : value(0)
{
  normalise();
}

Rational::Rational(const long _value) : value(_value)
{
  normalise();
}

Rational::Rational(const mp::Rational& _value) : value(_value)
{
  normalise();
}

Rational::Rational(const long num, const long denom) : 
  value(num, denom)
{
  assert(denom != 0);
  normalise();
}

std::size_t Rational::nops() const
{
  return 0;
}

void Rational::write(std::ostream& o) const
{
  o << value;
}

Expr Rational::derivative(const Symbol& s, Expr::optional_derivative_cache cache) const
{
  return zero();
}

bool Rational::depends(const std::set<Symbol>& symbols) const
{
  return false;
}

bool Rational::operator<(const Rational& n) const
{
  return value < n.value;
}

bool Rational::operator==(const long n) const
{
  return value == n;
}

bool Rational::operator==(const Rational& n) const
{
  return value == n.value;
}

Expr Rational::integrate(const Symbol& s, const unsigned flags) const
{
  return Sum::rational_multiple(s, *this).clone();
}

std::size_t Rational::untypedHash() const
{
  std::size_t hash = 0x161f15c2;
  cfd::util::hash_accum(hash, value);
  return hash;
}

Expr Rational::subs(const Expr::subst_map& map, const unsigned flags) const
{
  return clone();
}

Float Rational::eval(const Expr::subst_map& map) const
{
  return toFloat();
}

void Rational::accept(NumericExpressionVisitor<Symbol>& v) const
{
  v.visitConstant(getNumerator());

  if (getDenominator() != 1)
  {
    v.visitConstant(getDenominator());
    v.visitExponent(-1);
    v.postProduct(2);
  }
}

Rational Rational::operator-() const
{
  return Rational(-value);
}

Rational& Rational::operator+=(const Rational& r)
{
  value += r.value;
  normalise();
  return *this;
}

Rational& Rational::operator-=(const Rational& r)
{
  value -= r.value;
  normalise();
  return *this;
}

Rational& Rational::operator*=(const Rational& r)
{
  value *= r.value;
  normalise();
  return *this;
}

Rational& Rational::operator/=(const Rational& r)
{
  value /= r.value;
  normalise();
  return *this;
}

mp::Integer Rational::getNumerator() const
{
  return value.getNumerator();
}

mp::Integer Rational::getDenominator() const
{
  return value.getDenominator();
}

void Rational::normalise()
{
  invalidateHash();
}

Float Rational::toFloat() const
{
  return Float(value);
}

Rational Rational::reciprocal() const
{
  return value.reciprocal();
}

Rational Rational::gcd(const Rational& a, const Rational& b)
{
  using mp::Integer;

  if (a == 0)
    return b.abs();
  else if (b == 0)
    return a.abs();

  const Integer numerator = Integer::gcd(a.getNumerator(), b.getNumerator());
  const Integer denominator = Integer::gcd(a.getDenominator(), b.getDenominator());
  return Rational(mp::Rational(numerator, denominator));
}

Rational& Rational::operator++()
{
  ++value;
  normalise();
  return *this;
}

Rational& Rational::operator--()
{
  --value;
  normalise();
  return *this;
}

Rational Rational::abs() const
{
  return Rational(mp::abs(value));
}

Rational abs(const Rational& r)
{
  return r.abs();
}

Rational Rational::pow(const int exponent) const
{
  return Rational(mp::pow(value, exponent));
}

Rational pow(const Rational& r, const int exponent)
{
  return r.pow(exponent);
}

Expr Rational::extractMultiplier(Rational& coeff) const
{
  if (*this == 1)
  {
    return clone();
  }
  else
  {
    coeff *= *this;
    return Rational(1);
  }
}

Expr Rational::extractPolynomials(ExtractedExpressions& extracted) const
{
  return clone();
}

}

}
