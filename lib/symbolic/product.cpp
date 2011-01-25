#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/number.hpp>
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
  return Expr(new Number(1));
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
    ++newTerm[Number(d.second)];
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

    if (exponent == 0)
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
  return Number(0);
}

Expr Product::simplify() const
{
  const Number zero(0);
  const Expr simplified = PairSeq<Product>::simplify();
  const Basic& basic = simplified.internal();

  if (PairSeq<Product>::getType(basic) == PairSeq<Product>::getType(*this))
  {
    const Product& p = static_cast<const Product&>(basic);

    BOOST_FOREACH(const TermMap::value_type termPair, std::make_pair(p.begin(), p.end()))
    {
      // Eliminate 0^n where n>0.
      if (termPair.first == zero && termPair.second > 0)
        return Number(0);
    }
  }
  
  return simplified;
}

}

}
