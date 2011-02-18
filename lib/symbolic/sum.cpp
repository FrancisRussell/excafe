#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
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

    if (std::abs(current->second) != 1)
      o << std::abs(current->second) << "*";

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

  BOOST_FOREACH(const TermMap::value_type d, std::make_pair(begin(), end()))
  {
    d.first.accept(v);
    if (d.second != 1)
    {
      v.visitConstant(d.second);
      v.postProduct(2);
    }
  }

  v.postSummation(terms.size()+1);
}

Float Sum::eval(const Expr::subst_map& map) const
{
  Float result(getOverall().toFloat());

  BOOST_FOREACH(const TermMap::value_type d, std::make_pair(begin(), end()))
  {
    const Float evaluated = d.first.eval(map);
    result += evaluated * Float(d.second);;
  }
  return result;
}

void Sum::combineOverall(Rational& overall, const Rational& other)
{
  overall += other;
}

Rational Sum::applyCoefficient(const Rational& value, const int coefficient)
{
  return value * Rational(coefficient);
}


}

}
