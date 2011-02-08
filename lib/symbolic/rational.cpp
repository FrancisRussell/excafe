#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <boost/functional/hash.hpp>
#include <ostream>

namespace cfd
{

namespace symbolic
{

Rational::Rational(const int _value) : numerator(_value), denominator(1)
{
}

Rational::Rational(const int num, const int denom) : 
  numerator(num), denominator(denom)
{
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
  if (numerator < n.numerator)
    return true;
    
  if (numerator == n.numerator
      && denominator < n.denominator)
    return true;

  return false;
}

bool Rational::operator==(const Rational& n) const
{
  return numerator == n.numerator
         && denominator == n.denominator;
}

Expr Rational::integrate(const Symbol& s) const
{
  return Expr(new Product(*this, s));
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

Expr Rational::eval() const
{
  return Float(numerator)/Float(denominator);
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

}

}
