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

Sum Sum::operator+(const Expr& e) const
{
  TermMap map(begin(), end());
  ++map[e];
  return Sum(overall, map);
}

Expr Sum::derivative(const Symbol& s) const
{
  TermMap newTerms;
  BOOST_FOREACH(const TermMap::value_type& e, std::make_pair(begin(), end()))
  {
    newTerms.insert(std::make_pair(e.first.derivative(s), e.second));
  }
  return Sum(null(), newTerms);
}

Expr Sum::integrate(const Symbol& s) const
{
  TermMap dependentTerms;
  TermMap independentTerms;

  BOOST_FOREACH(const TermMap::value_type& e, std::make_pair(begin(), end()))
  {
    if (e.first.has(s))
      dependentTerms[e.first.integrate(s)] += e.second;
    else
      independentTerms[e.first] += e.second;
  }
  return Sum(null(), dependentTerms) + Product::mul(Sum(overall, independentTerms), s);
}

Sum Sum::expandedProduct(const Sum& other) const
{
  const Sum withoutOverallThis = this->withoutOverall();
  const Sum withoutOverallOther = other.withoutOverall();

  TermMap newTerms;
  BOOST_FOREACH(const TermMap::value_type& a, withoutOverallThis)
  {
    BOOST_FOREACH(const TermMap::value_type& b, withoutOverallOther)
    {
      newTerms[Product::mul(a.first, b.first).simplify()] += a.second*b.second;
    }
  }

  return Sum(null(), newTerms);
}

void Sum::accept(NumericExpressionVisitor<Symbol>& v) const
{
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

  v.postSummation(terms.size()+1);
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
  BOOST_FOREACH(const TermMap::value_type& d, terms)
  {
    result = Rational::gcd(result, d.second);
  }
  return result;
}

void Sum::extractMultipliers(Rational& overall, TermMap& map)
{
  TermMap newTermMap;
  BOOST_FOREACH(const TermMap::value_type& term, map)
  {
    Rational coefficient = term.second;
    const Expr e = term.first.internal().extractMultiplier(coefficient);
    newTermMap[e] += coefficient;
  }
  map.swap(newTermMap);
}


Expr Sum::extractMultiplier(Rational& coeff) const
{
  const Expr simplified = this->simplify();
  const Basic& basic = simplified.internal();

  if (is_a<Sum>(basic))
  {
    const Sum& sum = convert_to<Sum>(basic);
    const Rational multiplier = sum.findMultiplier();

    if (multiplier == Rational(1))
    {
      return simplified;
    }
    else
    {
      coeff *= multiplier;

      TermMap newTerms(sum.terms);
      BOOST_FOREACH(TermMap::value_type& d, newTerms)
      {
        d.second /= multiplier;
      }

      return Sum(sum.overall / multiplier, newTerms).clone();
    }
  }
  else
  {
    return basic.extractMultiplier(coeff);
  }
}

}

}
