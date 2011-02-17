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
  return Sum(map);
}

Expr Sum::derivative(const Symbol& s) const
{
  TermMap newTerms;
  BOOST_FOREACH(const TermMap::value_type& e, std::make_pair(begin(), end()))
  {
    newTerms.insert(std::make_pair(e.first.derivative(s), e.second));
  }
  return Sum(newTerms);
}

Expr Sum::integrate(const Symbol& s) const
{
  TermMap dependentTerms;
  TermMap independentTerms;
  BOOST_FOREACH(const TermMap::value_type& e, std::make_pair(begin(), end()))
  {
    if (e.first.has(s))
      dependentTerms.insert(std::make_pair(e.first.integrate(s), e.second));
    else
      independentTerms.insert(e);
  }
  return Sum(dependentTerms) + Product(Sum(independentTerms), s);
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
      newTerms[Product(a.first, b.first).simplify()] += a.second*b.second;
    }
  }

  return Sum(newTerms);
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

Expr Sum::eval() const
{
  Float floatPart(getOverall().toFloat());
  TermMap sum;

  BOOST_FOREACH(const TermMap::value_type d, std::make_pair(begin(), end()))
  {
    const Expr evaluated = d.first.eval();
    const Basic& value = evaluated.internal();
    if (getType(value) == getType(floatPart))
    {
      floatPart += Float(d.second) * static_cast<const Float&>(value);
    }
    else
    {
      sum[d.first] += d.second;
    }
  }
  
  if (sum.empty())
  {
    return floatPart;
  }
  else
  {
    ++sum[floatPart];
    return Sum(sum);
  }
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
