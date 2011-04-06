#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/exception.hpp>
#include <map>
#include <utility>
#include <cassert>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>

namespace cfd
{

namespace symbolic
{

Product Product::pow(const Expr& base, const int exponent)
{
  LazyTermMap terms;
  (*terms)[base] += exponent;
  return Product(null(), terms);
}

Product Product::div(const Expr& a, const Expr& b)
{
  LazyTermMap terms;
  ++(*terms)[a];
  --(*terms)[b];
  return Product(null(), terms);
}

Product Product::mul(const Expr& a, const Expr& b)
{
  LazyTermMap terms;
  ++(*terms)[a];
  ++(*terms)[b];
  return Product(null(), terms);
}

Rational Product::null()
{
  return Rational(1);
}

void Product::write(std::ostream& o) const
{
  if (begin() == end())
  {
    o << getOverall();
    return;
  }

  const_iterator current(begin());

  o << "(";

  if (getOverall() != null())
    o << getOverall() << "*";

  while(current != end())
  {
    const bool expIsOne = current->second == 1;
    o << (expIsOne ? "" : "(") << current->first << (expIsOne ? "" : ")");

    if (!expIsOne) 
      o << "^" << "{" << current->second << "}";

    ++current;

    if (current != end())
      o << "*";
  }
  o << ")";
}

Expr Product::derivative(const Symbol& s) const
{
  Sum summation;

  BOOST_FOREACH(const TermMap::value_type& d, *this)
  {
    LazyTermMap newTerm(getTerms());
    newTerm->erase(d.first);
    ++(*newTerm)[Rational(d.second)];
    ++(*newTerm)[d.first.derivative(s)];
    (*newTerm)[d.first]+=d.second-1;
    summation += Product(getOverall(), newTerm);
  }
  
  return summation;
}

Expr Product::integrate_internal(const Symbol& s) const
{
  LazyTermMap independent;
  TermMap dependent;

  /* We factor into products dependent and not dependent on s */
  BOOST_FOREACH(const TermMap::value_type& d, *this)
  {
    if (d.first.depends(s))
      dependent.insert(d);
    else
      independent->insert(d);
  }

  Expr dependentIntegral;
  if (dependent.empty())
  {
    /* Integration of 1 */
    dependentIntegral = s;
  }
  else if (dependent.size() == 1)
  {
    const Expr expr = dependent.begin()->first;
    const int exponent = dependent.begin()->second;

    if (exponent < 0)
      CFD_EXCEPTION("Cannot integrate functions involving variable raised to negative exponents.");

    if (expr == s)
    {
      /* Directly handle variables raised to an exponent */
      dependentIntegral = Product::mul(Rational(1, exponent+1), Product::pow(s, exponent+1));
    }
    else if (exponent == 0)
    {
      /* Integral of f^0 (if we've failed to eliminate them for some reason)*/
      dependentIntegral = s;
    }
    else if (exponent == 1)
    {
      /* Integration of f^1*/
      dependentIntegral = expr.integrate_internal(s); 
    }
    else
    {
      /* Integration of f^n where n>1 */
      const int exp1 = exponent / 2;
      const int exp2 = exponent - exp1;
      dependentIntegral = integrate(pow(expr, exp1), pow(expr, exp2), s);
    }
  }
  else
  {
    const TermMap::iterator pivot = 
      boost::next(dependent.begin(), dependent.size()/2);

    const TermMap first(dependent.begin(), pivot);
    const TermMap second(pivot, dependent.end());

    dependentIntegral = integrate(Product(null(), first), Product(null(), second), s);
  }

  return Product(getOverall(), independent) * dependentIntegral;
}

Expr Product::integrate(const Expr& a, const Expr& b, const Symbol& s)
{
  const Rational zero(0);
  int sign = 1;
  Sum result;

  Expr u = a;
  Expr v = b.integrate_internal(s);

  while (u != zero)
  {
    result += Sum::rational_multiple(Product::mul(u, v), Rational(sign));
    u = u.derivative(s);
    v = v.integrate_internal(s);
    sign *= -1;
  }

  return result.simplify();
}

void Product::accept(NumericExpressionVisitor<Symbol>& v) const
{
  const bool hasOverall = (getOverall() != null());
  if (hasOverall)
    getOverall().accept(v);

  BOOST_FOREACH(const TermMap::value_type& d, *this)
  {
    d.first.accept(v);
    if (d.second != 1)
      v.visitExponent(d.second);
  }

  v.postProduct(getTerms().size() + (hasOverall ? 1 : 0));
}

Float Product::eval(const Expr::subst_map& map) const
{
  Float result(getOverall().toFloat());

  BOOST_FOREACH(const TermMap::value_type& d, *this)
  {
    const Float evaluated = d.first.eval(map);
    result *= symbolic::pow(evaluated, d.second);
  }

  return result;
}

void Product::combineOverall(Rational& overall, const Rational& other)
{
  overall *= other;
}

Rational Product::applyCoefficient(const Rational& value, const int coefficient)
{
  return symbolic::pow(value, coefficient);
}

Product Product::extractMultipliers() const
{
  Rational overall = getOverall();
  LazyTermMap newTermMap;

  BOOST_FOREACH(const TermMap::value_type& term, getTerms())
  {
    Rational multiplier = null();
    const Expr newTerm = term.first.internal().extractMultiplier(multiplier);
    (*newTermMap)[newTerm] += term.second;
    overall *= symbolic::pow(multiplier, term.second);
  }

  // Remove all terms if overall coefficient is zero.
  if (overall == 0)
    newTermMap->clear();
  
  return Product(overall, newTermMap);
}

Product& Product::operator*=(const Product& p)
{
  this->combine(p);
  return *this;
}

Expr Product::extractMultiplier(Rational& coeff) const
{
  if (this->getRewriteState() == NORMALISED_AND_EXTRACTED)
    return clone();

  const Product product  = (this->getRewriteState() == NORMALISED ? *this : this->getNormalised());
  coeff *= product.overall;
  return constructSimplifiedExpr(null(), product.terms, NORMALISED_AND_EXTRACTED);
}

}

}
