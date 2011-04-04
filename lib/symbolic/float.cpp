#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <boost/functional/hash.hpp>
#include <cmath>

namespace cfd
{

namespace symbolic
{

Float Float::fromFraction(const long numerator, const long denominator)
{
  return Float(static_cast<double>(numerator)/denominator);
}

Float::Float(const double _value) : value(_value)
{
}

std::size_t Float::nops() const
{
  return 0;
}

void Float::write(std::ostream& o) const
{
  o << value;
}

Expr Float::derivative(const Symbol& s) const
{
  return Expr(new Rational(0));
}

bool Float::has(const Expr& e) const
{
  return e == *this;
}

bool Float::operator<(const Float& n) const
{
  return value < n.value;
}

bool Float::operator==(const Float& n) const
{
  return value == n.value;
}

Float& Float::operator+=(const Float& n)
{
  invalidateHash();
  value += n.value;
  return *this;
}

Float& Float::operator-=(const Float& n)
{
  invalidateHash();
  value -= n.value;
  return *this;
}

Float& Float::operator/=(const Float& n)
{
  invalidateHash();
  value /= n.value;
  return *this;
}

Float& Float::operator*=(const Float& n)
{
  invalidateHash();
  value *= n.value;
  return *this;
}

Expr Float::integrate_internal(const Symbol& s) const
{
  return Product::mul(*this, s);
}

std::size_t Float::untypedHash() const
{
  std::size_t result = 0x2c6831da;
  boost::hash_combine(result, value);
  return result;
}

Expr Float::subs(const Expr::subst_map& map) const
{
  return clone();
}

Float Float::eval(const Expr::subst_map& map) const
{
  return *this;
}

Expr Float::simplify() const
{
  const long multiplier = 2 << 6;
  const long truncated = static_cast<long>(multiplier * value);
  if (multiplier * value == truncated)
    return Rational(truncated, multiplier);
  else
    return clone();
}

void Float::accept(NumericExpressionVisitor<Symbol>& v) const
{
  v.visitConstant(value);
}

double Float::toDouble() const
{
  return value;
}

Float Float::pow(const int exponent) const
{
  return Float(std::pow(value, exponent));
}

Float pow(const Float& f, const int exponent)
{
  return f.pow(exponent);
}

}

}
