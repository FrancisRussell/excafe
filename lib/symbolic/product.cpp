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
    TermMap newTerm(begin(), end());
    newTerm.erase(d.first);
    ++newTerm[Rational(d.second)];
    ++newTerm[d.first.derivative(s)];
    newTerm[d.first]+=d.second-1;
    summation = summation + Product(getOverall(), newTerm);
  }
  
  return summation;
}

Expr Product::integrate(const Symbol& s) const
{
  TermMap independent;
  TermMap dependent;

  /* We factor into products dependent and not dependent on s */
  BOOST_FOREACH(const TermMap::value_type& d, *this)
  {
    if (d.first.has(s))
      dependent.insert(d);
    else
      independent.insert(d);
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
      dependentIntegral = expr.integrate(s); 
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
  Expr result = zero;

  Expr u = a;
  Expr v = b.integrate(s);

  while (u != zero)
  {
    result = result + Sum::rational_multiple(Product::mul(u, v), Rational(sign));
    u = u.derivative(s);
    v = v.integrate(s);
    sign *= -1;
  }

  return result.simplify();
}

Expr Product::simplify() const
{
  const Rational zero(0);
  const Expr simplified = PairSeq<Product, int>::simplify();
  const Basic& basic = simplified.internal();

  if (is_a<Product>(basic))
  {
    const Product& p = convert_to<Product>(basic);

    if (p.getOverall() == zero)
      return zero;
  }
  
  return simplified;
}

void Product::accept(NumericExpressionVisitor<Symbol>& v) const
{
  getOverall().accept(v);
  BOOST_FOREACH(const TermMap::value_type& d, *this)
  {
    d.first.accept(v);
    if (d.second != 1)
      v.visitExponent(d.second);
  }

  v.postProduct(terms.size()+1);
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

void Product::extractMultipliers(Rational& overall, TermMap& map)
{
  TermMap newTermMap;
  BOOST_FOREACH(const TermMap::value_type& term, map)
  {
    Rational multiplier = null();
    const Expr newTerm = term.first.internal().extractMultiplier(multiplier);
    newTermMap[newTerm] += term.second;
    overall *= symbolic::pow(multiplier, term.second);
  }
  map.swap(newTermMap);
}

Product& Product::operator*=(const Product& p)
{
  this->combine(p);
  return *this;
}

Expr Product::extractMultiplier(Rational& coeff) const
{
  const Expr simplified = this->simplify();
  const Basic& basic = simplified.internal();

  if (is_a<Product>(basic))
  {
    const Product& product = convert_to<Product>(basic);
    if (product.getOverall() == null())
    {
      return simplified;
    }
    else
    {
      coeff *= product.overall;
      return Product(null(), product.terms).clone();
    }
  }
  else
  {
    return basic.extractMultiplier(coeff);
  }
}

}

}
