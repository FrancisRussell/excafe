#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/flags.hpp>
#include <simple_cfd/symbolic/collect_visitor.hpp>
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

Product Product::constant(const Rational& r)
{
  LazyTermMap terms;
  return Product(r, terms);
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

/*
  Generalised power rule: (f^g)' (e^(g*ln f))' = f^g * (f'gf^-1 + g'ln f)
  Product rule (for 3 functions): (fgh)' = f'gh + fg'h + fgh'
*/

Expr Product::derivative(const Symbol& s, Expr::optional_derivative_cache cache) const
{
  Sum summation;

  BOOST_FOREACH(const TermMap::value_type& d, *this)
  {
    const Expr termDerivative = d.first.derivative(s, cache);

    if (termDerivative != Rational::zero())
    {
      const Rational termOverall = getOverall() * d.second;
      LazyTermMap newTerm(getTerms());

      if (d.second == 1)
        newTerm->erase(d.first);
      else
        --(*newTerm)[d.first];

      ++(*newTerm)[termDerivative];
      summation += constructSimplifiedExpr(termOverall, newTerm, NON_NORMALISED);
    }
  }
  
  return summation;
}

Expr Product::integrateComplex(const LazyTermMap& terms, const Symbol& s, const unsigned flags)
{
  if ((flags & Flags::DO_NOT_COLLECT) == 0)
  {
    const Expr product = Product(null(), terms).clone();

    CollectVisitor collectVisitor(s);
    product.accept(collectVisitor);
    const Expr collected = collectVisitor.getResult();

    return collected.integrate(s, flags | Flags::DO_NOT_COLLECT);
  }
  else if (terms->size() == 1)
  {
    const Expr expr = terms->begin()->first;
    const int exponent = terms->begin()->second;

    const int exp1 = exponent / 2;
    const int exp2 = exponent - exp1;

    return integrate(pow(expr, exp1), pow(expr, exp2), s, flags);
  }
  else
  {
    const TermMap::const_iterator pivot = 
      boost::next(terms->begin(), terms->size()/2);

    LazyTermMap first, second;
    first->insert(terms->begin(), pivot);
    second->insert(pivot, terms->end());

    return integrate(Product(null(), first), Product(null(), second), s, flags);
  }
}

Expr Product::integrate(const Symbol& s, const unsigned flags) const
{
  LazyTermMap dependent;
  LazyTermMap independent;

  std::set<Symbol> symbols;
  symbols.insert(s);
  this->extractDependent(symbols, *dependent, *independent);

  Expr dependentIntegral;
  if (dependent->empty())
  {
    /* Integration of 1 */
    dependentIntegral = s;
  }
  else if (dependent->size() == 1)
  {
    const Expr expr = dependent->begin()->first;
    const int exponent = dependent->begin()->second;

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
      dependentIntegral = expr.integrate(s, flags); 
    }
    else
    {
      /* Integration of f^n where n>1 */
      dependentIntegral = integrateComplex(dependent, s, flags);
    }
  }
  else
  {
    dependentIntegral = integrateComplex(dependent, s, flags);
  }

  ++(*independent)[dependentIntegral];
  return constructSimplifiedExpr(getOverall(), independent, NON_NORMALISED);
}

Expr Product::integrate(const Product& a, const Product& b, const Symbol& s, const unsigned flags)
{
  int sign = 1;
  Sum result;

  Expr u = a;
  Expr v = b.integrate(s, flags);

  while (u != Rational::zero())
  {
    result += Sum::rational_multiple(Product::mul(u, v), Rational(sign));
    u = u.derivative(s);
    v = v.integrate(s, flags);
    sign *= -1;
  }

  return result;
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

Expr Product::integrate(const Expr::region_t& region, const unsigned flags) const
{
  LazyTermMap dependent;
  LazyTermMap independent;

  const std::set<Symbol> symbols = region.getVariables();
  this->extractDependent(symbols, *dependent, *independent);

  Expr integrated = Product(null(), dependent).clone();
  if ((flags & Flags::DO_NOT_COLLECT) == 0)
  {
    CollectVisitor collectVisitor(region.getVariables());
    integrated.accept(collectVisitor);
    integrated = collectVisitor.getIntegratedResult(region, flags);
  }
  else
  {
    BOOST_FOREACH(const Expr::region_t::value_type& interval, region)
    {
      const Symbol& variable = interval.first;
      integrated = integrated.integrate(variable, flags);
  
      Expr::subst_map lower, upper;
      lower[variable] = interval.second.first;
      upper[variable] = interval.second.second;
      integrated = (integrated.subs(upper) - integrated.subs(lower)).simplify();
    }
  }

  ++(*independent)[integrated];
  return constructSimplifiedExpr(getOverall(), independent, NON_NORMALISED);
}

}

}
