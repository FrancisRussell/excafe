#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/make_expr_from.hpp>
#include <simple_cfd/numeric/cast.hpp>
#include <simple_cfd/util/hash.hpp>
#include <cln/real.h>
#include <set>
#include <cmath>

namespace cfd
{

namespace symbolic
{

Float::Float() : value(0.0)
{
}

Float::Float(const double _value) : value(_value)
{
}

Float::Float(const mp::Rational& _value) : 
  value(_value.toDouble())
{
}

Float::Float(const cln::cl_R& _value) : 
  value(cln::double_approx(_value))
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

Expr Float::derivative(const Symbol& s, Expr::optional_derivative_cache cache) const
{
  return Rational::zero();
}

bool Float::depends(const std::set<Symbol>& symbols) const
{
  return false;
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

Expr Float::integrate(const Symbol& s, const unsigned flags) const
{
  return Product::mul(*this, s);
}

std::size_t Float::untypedHash() const
{
  std::size_t result = 0x2c6831da;
  cfd::util::hash_accum(result, value);
  return result;
}

Expr Float::subs(const Expr::subst_map& map, const unsigned flags) const
{
  return clone();
}

Float Float::eval(const Expr::subst_map& map) const
{
  return *this;
}

Expr Float::simplify() const
{
  Rational rational;
  if (asRational(rational))
    return rational;
  else
    return clone();
}

bool Float::asRational(Rational& r) const
{
  const long multiplier = 1 << 8;
  const long truncated = static_cast<long>(multiplier * value);
  if (multiplier * value == truncated)
  {
    r = Rational(truncated, multiplier);
    return true;
  }
  else
  {
    return false;
  }
}

Expr Float::extractMultiplier(Rational& coeff) const
{
  Rational rational;
  if (asRational(rational))
  {
    coeff *= rational;
    return Rational::one();
  }
  else
  {
    return clone();
  }
}

void Float::accept(NumericExpressionVisitor<Symbol>& v) const
{
  v.visitConstant(cfd::numeric_cast<NumericExpressionVisitor<Symbol>::float_t>(value));
}

double Float::toDouble() const
{
  return value;
}

Float Float::pow(const int exponent) const
{
  return Float(std::pow(value, exponent));
}

Float Float::abs() const
{
  return Float(std::abs(value));
}

Float pow(const Float& f, const int exponent)
{
  return f.pow(exponent);
}

Expr Float::extractPolynomials(ExtractedExpressions& extracted) const
{
  return clone();
}

}

}
