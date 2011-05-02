#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <map>
#include <utility>
#include <boost/foreach.hpp>

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

Sum Sum::constant(const Rational& r)
{
  LazyTermMap terms;
  return Sum(r, terms);
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

Sum& Sum::operator/=(const Rational& r)
{
  if (r == 1)
    return *this;

  this->invalidateHash();
  this->overall /= r;

  BOOST_FOREACH(TermMap::value_type& e, getTerms())
    e.second /= r;

  return *this;
}

Expr Sum::derivative(const Symbol& s) const
{
  LazyTermMap newTerms;

  BOOST_FOREACH(const TermMap::value_type& e, getTerms())
  {
    const Expr d = e.first.derivative(s);

    if (d != Rational::zero())
      (*newTerms)[d] += e.second;
  }

  return constructSimplifiedExpr(Rational(0), newTerms, NON_NORMALISED);
}

Expr Sum::integrate(const Symbol& s, const unsigned flags) const
{
  LazyTermMap dependentTerms;
  LazyTermMap independentTerms;

  BOOST_FOREACH(const TermMap::value_type& e, getTerms())
  {
    if (e.first.depends(s))
      (*dependentTerms)[e.first.integrate(s, flags)] += e.second;
    else
      (*independentTerms)[e.first] += e.second;
  }

  Sum result(null(), dependentTerms);
  result += Product::mul(Sum(overall, independentTerms), s);
  return result;
}

Sum Sum::expandedProduct(const Sum& other) const
{
  const Sum withoutOverallThis = this->withoutOverall();
  const Sum withoutOverallOther = other.withoutOverall();

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
  std::size_t negativeCount = 0;
  if (this->getOverall() < 0) 
    ++negativeCount;

  Rational result = this->getOverall();
  BOOST_FOREACH(const TermMap::value_type& d, getTerms())
  {
    if (d.second < 0)
      ++negativeCount;

    result = Rational::gcd(result, d.second);
  }

  // If more than half of the terms in the sum have a negative
  // co-efficient, we negate the extracted multiplier. This leads to a
  // consistent normal form except when there are the same number of
  // positive and negative terms.
  const std::size_t numTerms = getTerms().size() + (getOverall() == 0 ? 0 : 1);
  if (negativeCount > numTerms/2)
    result = -result;

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
  if (this->getRewriteState() == NORMALISED_AND_EXTRACTED)
    return clone();

  const Sum sum = (this->getRewriteState() == NORMALISED ? *this : this->getNormalised());
  const Rational multiplier = sum.findMultiplier();

  coeff *= multiplier;

  // Even though gcd(0,n) == |n|, we still need to handle the all-zero case.
  if (multiplier == 0)
    return Rational::one();

  LazyTermMap newTerms = sum.terms;

  // We can avoid forcing a copy of the TermMap if multiplier==1.
  if (multiplier != 1)
  {
    BOOST_FOREACH(TermMap::value_type& d, *newTerms)
    {
      d.second /= multiplier;
    }
  }

  return constructSimplifiedExpr(sum.overall / multiplier, newTerms, NORMALISED_AND_EXTRACTED);
}

Expr Sum::integrate(const Expr::region_t& region, const unsigned flags) const
{
  LazyTermMap resultTerms;
  BOOST_FOREACH(const TermMap::value_type& e, getTerms())
  {
    (*resultTerms)[e.first.integrate(region, flags)] += e.second;
  }

  return constructSimplifiedExpr(getOverall() * region.getVolume(), resultTerms, NON_NORMALISED);
}

}

}
