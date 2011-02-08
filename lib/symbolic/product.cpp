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

Expr Product::null() const
{
  return Expr(new Rational(1));
}

void Product::write(std::ostream& o) const
{
  if (begin() == end())
  {
    o << "1.0";
    return;
  }

  const_iterator current(begin());

  o << "(";
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

  BOOST_FOREACH(const TermMap::value_type d, std::make_pair(begin(), end()))
  {
    TermMap newTerm(begin(), end());
    newTerm.erase(d.first);
    ++newTerm[Rational(d.second)];
    ++newTerm[d.first.derivative(s)];
    newTerm[d.first]+=d.second-1;
    summation = summation + Product(newTerm);
  }
  
  return summation;
}

Expr Product::integrate(const Symbol& s) const
{
  TermMap independent;
  TermMap dependent;

  /* We factor into products dependent and not dependent on s */
  BOOST_FOREACH(const TermMap::value_type d, std::make_pair(begin(), end()))
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
      dependentIntegral = Product(Rational(1, exponent+1), Product::pow(s, exponent+1));
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

    dependentIntegral = integrate(Product(first), Product(second), s);
  }

  return Product(independent) * dependentIntegral;
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
    result = result + Sum::integer_multiple(Product(u, v), sign);
    u = u.derivative(s);
    v = v.integrate(s);
    sign *= -1;
  }

  return result.simplify();
}

Expr Product::simplify() const
{
  const Rational zero(0);
  const Expr simplified = PairSeq<Product>::simplify();
  const Basic& basic = simplified.internal();

  if (PairSeq<Product>::getType(basic) == PairSeq<Product>::getType(*this))
  {
    const Product& p = static_cast<const Product&>(basic);

    BOOST_FOREACH(const TermMap::value_type termPair, std::make_pair(p.begin(), p.end()))
    {
      // Eliminate 0^n where n>0.
      if (termPair.first == zero && termPair.second > 0)
        return Rational(0);
    }
  }
  
  return simplified;
}

void Product::accept(NumericExpressionVisitor<Symbol>& v) const
{
  BOOST_FOREACH(const TermMap::value_type d, std::make_pair(begin(), end()))
  {
    d.first.accept(v);
    if (d.second != 1)
      v.visitExponent(d.second);
  }

  v.postProduct(terms.size());
}

Expr Product::eval() const
{
  Float floatPart(1.0);
  TermMap prod;

  BOOST_FOREACH(const TermMap::value_type d, std::make_pair(begin(), end()))
  {
    const Expr evaluated = d.first.eval();
    const Basic& value = evaluated.internal();
    if (getType(value) == getType(floatPart))
    {
      if (d.second > 0)
      {
        for (int i=0; i<d.second; ++i)
          floatPart *= static_cast<const Float&>(value);
      }
      else
      {
        for (int i=0; i<-d.second; ++i)
          floatPart /= static_cast<const Float&>(value);
      }
    }
    else
    {
      prod[d.first] += d.second;
    }
  }
  
  if (prod.empty())
  {
    return floatPart;
  }
  else
  {
    ++prod[floatPart];
    return Product(prod);
  }
}

}

}
