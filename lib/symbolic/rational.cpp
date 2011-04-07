#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/util/hash.hpp>
#include <ostream>
#include <cassert>
#include <cmath>
#include <set>

namespace cfd
{

namespace symbolic
{

Rational::Rational() : numerator(0), denominator(1)
{
  normalise();
}

Rational::Rational(const value_type _value) : numerator(_value), denominator(1)
{
  normalise();
}

Rational::Rational(const value_type num, const value_type denom) : 
  numerator(num), denominator(denom)
{
  assert(denominator != 0);
  normalise();
}

std::size_t Rational::nops() const
{
  return 0;
}

void Rational::write(std::ostream& o) const
{
  const bool isInteger = (denominator == 1);

  if (!isInteger)
    o << "(";

  o << numerator;

  if (!isInteger)
    o << "/" << denominator << ")";
}

Expr Rational::derivative(const Symbol& s) const
{
  return Expr(new Rational(0));
}

bool Rational::depends(const std::set<Symbol>& symbols) const
{
  return false;
}

bool Rational::operator<(const Rational& n) const
{
  if (*this == n)
    return false;
  else
    return numerator*n.denominator 
           < n.numerator*denominator;
}

bool Rational::operator==(const Rational& n) const
{
  return numerator == n.numerator
         && denominator == n.denominator;
}

Expr Rational::integrate_internal(const Symbol& s) const
{
  return Sum::rational_multiple(s, *this).clone();
}

std::size_t Rational::untypedHash() const
{
  std::size_t hash = 0x161f15c2;
  cfd::util::hash_accum(hash, numerator);
  cfd::util::hash_accum(hash, denominator);
  return hash;
}

Expr Rational::subs(const Expr::subst_map& map) const
{
  return clone();
}

Float Rational::eval(const Expr::subst_map& map) const
{
  return toFloat();
}

void Rational::accept(NumericExpressionVisitor<Symbol>& v) const
{
  v.visitConstant(numerator);

  if (denominator != 1)
  {
    v.visitConstant(denominator);
    v.visitExponent(-1);
    v.postProduct(2);
  }
}

Rational Rational::operator-() const
{
  return Rational(-numerator, denominator);
}

Rational& Rational::operator+=(const Rational& r)
{
  const value_type lcd = lcm(denominator, r.denominator);
  numerator = numerator * (lcd/denominator) + r.numerator * (lcd/r.denominator);
  denominator = lcd;
  normalise();
  return *this;
}

Rational& Rational::operator-=(const Rational& r)
{
  *this += -r;
  return *this;
}

Rational& Rational::operator*=(const Rational& r)
{
  const value_type gcd1 = gcd(numerator, r.denominator);
  const value_type gcd2 = gcd(r.numerator, denominator);

  numerator = (numerator/gcd1)*(r.numerator/gcd2);
  denominator = (denominator/gcd2)*(r.denominator/gcd1);
  normalise();
  return *this;
}

Rational& Rational::operator/=(const Rational& r)
{
  *this *= r.reciprocal();
  return *this;
}

Rational::value_type Rational::getNumerator() const
{
  return numerator;
}

Rational::value_type Rational::getDenominator() const
{
  return denominator;
}

void Rational::normalise()
{
  invalidateHash();

  if (numerator == 0)
    denominator = 1;

  // Short circuit for efficiency.
  if (denominator == 1)
    return;

  if (denominator < 0)
  {
    numerator = -numerator;
    denominator = -denominator;
  }

  const value_type factor = gcd(numerator, denominator);
  numerator /= factor;
  denominator /= factor;

  assert(denominator != 0);
}

Float Rational::toFloat() const
{
  return Float::fromFraction(numerator, denominator);
}

Rational Rational::reciprocal() const
{
  return Rational(denominator, numerator);
}

Rational::value_type Rational::gcd(value_type u, value_type v)
{
  u = std::abs(u);
  v = std::abs(v);

  if (u == 0 || v == 0)
    return u | v;

  unsigned shift;
  for(shift=0; ((u | v) & 1) == 0; ++shift)
  {
    u >>= 1;
    v >>= 1;
  }

  while ((u & 1) == 0)
    u >>= 1;

  do
  {
    while ((v & 1) == 0)
      v >>= 1;

    if (u > v)
    {
      const value_type tmp = u;
      u = v;
      v = tmp;
    }

    v -= u;
    v >>= 1;
  }
  while (v != 0);

  return u << shift;
}

Rational::value_type Rational::lcm(const value_type a, const value_type b)
{
  // Performing the division first reduces the size of the intermediate
  // value when gcd(a,b) > 1.

  return b*(a/gcd(a, b));
}

Rational Rational::gcd(const Rational& a, const Rational& b)
{
  if (a.numerator == 0)
    return abs(b);
  else if (b.numerator == 0)
    return abs(a);

  const value_type numerator = gcd(a.numerator, b.numerator);
  const value_type denominator = gcd(a.denominator, b.denominator);
  return Rational(numerator, denominator);
}

Rational& Rational::operator++()
{
  numerator+=denominator;
  normalise();
  return *this;
}

Rational& Rational::operator--()
{
  numerator-=denominator;
  normalise();
  return *this;
}

Rational abs(const Rational& r)
{
  if (r.getNumerator() < 0)
    return -r;
  else
    return r;
}

Rational pow(const Rational& r, const int exponent)
{
  const Rational multiplier = exponent >=0 ? r : r.reciprocal();
  const int repetitions = std::abs(exponent);

  Rational result(1);
  for(int i=0; i<repetitions; ++i)
    result *= multiplier;

  return result;
}

Expr Rational::extractMultiplier(Rational& coeff) const
{
  if (*this == Rational(1))
  {
    return clone();
  }
  else
  {
    coeff *= *this;
    return Rational(1);
  }
}

}

}
