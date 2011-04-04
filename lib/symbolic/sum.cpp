#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <map>
#include <utility>
#include <boost/foreach.hpp>
#include <boost/function.hpp>

namespace cfd
{

namespace symbolic
{

Sum Sum::sub(const Expr& a, const Expr& b)
{
  LazyTermMap terms;
  ++(*terms)[a];
  --(*terms)[b];
  return Sum(null(), terms);
}

Sum Sum::rational_multiple(const Expr& e, const Rational& n)
{
  LazyTermMap terms;
  (*terms)[e]+=n;
  return Sum(null(), terms);
}

Sum Sum::add(const Expr& a, const Expr& b)
{
  LazyTermMap terms;
  ++(*terms)[a];
  ++(*terms)[b];
  return Sum(null(), terms);
}

bool Sum::AlwaysTrue(const Expr& e)
{
  return true;
}

Rational Sum::null()
{
  return Rational(0);
}

void Sum::write(std::ostream& o) const
{
  if (begin() == end())
  {
    o << getOverall();
    return;
  }

  const_iterator current(begin());

  o << "(";

  if (getOverall() != null())
    o << getOverall() << "+";

  while(current != end())
  {
    o << (current->second < 0 ? "-" : "");

    if (abs(current->second) != 1)
      o << abs(current->second) << "*";

    o << current->first;
    ++current;

    if (current != end())
      o << " + ";
  }
  o << ")";
}

Sum& Sum::operator+=(const Expr& e)
{
  this->invalidateHash();
  ++getTerms()[e];
  return *this;
}

Expr Sum::derivative(const Symbol& s) const
{
  Sum result;
  BOOST_FOREACH(const TermMap::value_type& e, getTerms())
  {
    result.getTerms()[e.first.derivative(s)] += e.second;
  }
  return result;
}

Expr Sum::integrate_internal(const Symbol& s) const
{
  LazyTermMap dependentTerms;
  LazyTermMap independentTerms;

  BOOST_FOREACH(const TermMap::value_type& e, getTerms())
  {
    if (e.first.has(s))
      (*dependentTerms)[e.first.integrate_internal(s)] += e.second;
    else
      (*independentTerms)[e.first] += e.second;
  }

  Sum result(null(), dependentTerms);
  result += Product::mul(Sum(overall, independentTerms), s);
  return result;
}

Sum Sum::groupNonMatching(const Sum& sum, const boost::function<bool (const Expr&)>& predicate)
{
  if (predicate != AlwaysTrue)
  {
    const Sum withoutOverall = sum.withoutOverall();
    LazyTermMap matching, nonMatching;
  
    BOOST_FOREACH(const TermMap::value_type& term, withoutOverall)
    {
      (predicate(term.first) ? *matching : *nonMatching)[term.first] += term.second;
    }
  
    if (!nonMatching->empty())
      ++(*matching)[Sum(null(), nonMatching).clone()];
  
    return Sum(null(), matching);
  }
  else
  {
    return sum;
  }
}

Sum Sum::expandedProduct(const Sum& other, const boost::function<bool (const Expr&)>& predicate) const
{
  const Sum withoutOverallThis = groupNonMatching(*this, predicate).withoutOverall();
  const Sum withoutOverallOther = groupNonMatching(other, predicate).withoutOverall();

  LazyTermMap newTerms;
  BOOST_FOREACH(const TermMap::value_type& a, withoutOverallThis)
  {
    BOOST_FOREACH(const TermMap::value_type& b, withoutOverallOther)
    {
      (*newTerms)[Product::mul(a.first, b.first).simplify()] += a.second*b.second;
    }
  }

  return Sum(null(), newTerms);
}

void Sum::accept(NumericExpressionVisitor<Symbol>& v) const
{
  const bool hasOverall = (getOverall() != null());
  if (hasOverall)
    getOverall().accept(v);

  BOOST_FOREACH(const TermMap::value_type& d, *this)
  {
    d.first.accept(v);
    if (d.second != 1)
    {
      d.second.accept(v);
      v.postProduct(2);
    }
  }

  v.postSummation(getTerms().size() + (hasOverall ? 1 : 0));
}

Float Sum::eval(const Expr::subst_map& map) const
{
  Float result(getOverall().toFloat());

  BOOST_FOREACH(const TermMap::value_type& d, *this)
  {
    const Float evaluated = d.first.eval(map);
    result += evaluated * d.second.toFloat();
  }
  return result;
}

void Sum::combineOverall(Rational& overall, const Rational& other)
{
  overall += other;
}

Rational Sum::applyCoefficient(const Rational& value, const Rational& coefficient)
{
  return value * Rational(coefficient);
}

Rational Sum::findMultiplier() const
{
  Rational result = this->getOverall();
  BOOST_FOREACH(const TermMap::value_type& d, getTerms())
  {
    result = Rational::gcd(result, d.second);
  }

  // Even though gcd(0,n) == |n|, we still need to handle the all-zero case.
  if (result == 0)
    return Rational(1);
  else
    return result;
}

Sum Sum::extractMultipliers() const
{
  LazyTermMap newTermMap;

  BOOST_FOREACH(const TermMap::value_type& term, getTerms())
  {
    Rational coefficient = term.second;
    const Expr e = term.first.internal().extractMultiplier(coefficient);
    (*newTermMap)[e] += coefficient;
  }

  return Sum(getOverall(), newTermMap);
}

Sum& Sum::operator+=(const Sum& s)
{
  this->combine(s);
  return *this;
}

Expr Sum::extractMultiplier(Rational& coeff) const
{
  const Sum sum = this->getNormalised();
  const Rational multiplier = sum.findMultiplier();

  coeff *= multiplier;

  LazyTermMap newTerms(sum.getTerms());
  BOOST_FOREACH(TermMap::value_type& d, *newTerms)
  {
    d.second /= multiplier;
  }

  return constructSimplifiedExpr(sum.overall / multiplier, newTerms);
}

}

}
