#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <boost/functional/hash.hpp>
#include <ostream>
#include <cassert>
#include <cmath>

namespace cfd
{

namespace symbolic
{

Rational::Rational() : numerator(0), denominator(1)
{
  normalise();
}

Rational::Rational(const long _value) : numerator(_value), denominator(1)
{
  normalise();
}

Rational::Rational(const long num, const long denom) : 
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

bool Rational::has(const Expr& e) const
{
  return e == *this;
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

Expr Rational::integrate(const Symbol& s) const
{
  return Product::mul(*this, s).clone();
}

std::size_t Rational::untypedHash() const
{
  std::size_t hash = 0;
  boost::hash_combine(hash, numerator);
  boost::hash_combine(hash, denominator);
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
  numerator = numerator * r.denominator + r.numerator * denominator;
  denominator *= r.denominator;
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
  numerator *= r.numerator;
  denominator *= r.denominator;
  normalise();
  return *this;
}


Rational& Rational::operator/=(const Rational& r)
{
  numerator *= r.denominator;
  denominator *= r.numerator;
  normalise();
  return *this;
}

long Rational::getNumerator() const
{
  return numerator;
}

long Rational::getDenominator() const
{
  return denominator;
}

void Rational::normalise()
{
  invalidateHash();

  if (numerator == 0)
    denominator = 1;

  if (denominator < 0)
  {
    numerator = -numerator;
    denominator = -denominator;
  }

  const long factor = gcd(std::abs(numerator), std::abs(denominator));
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

long Rational::gcd(long u, long v)
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
      const unsigned tmp = u;
      u = v;
      v = tmp;
    }

    v -= u;
    v >>= 1;
  }
  while (v != 0);

  return u << shift;
}

Rational Rational::gcd(const Rational& a, const Rational& b)
{
  const long numerator = gcd(a.numerator, b.numerator);
  const long denominator = gcd(a.denominator, b.denominator);
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

}

}
