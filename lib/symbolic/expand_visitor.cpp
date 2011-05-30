#include <utility>
#include <cassert>
#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <simple_cfd/symbolic/expand_visitor.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/group.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/basic.hpp>

namespace cfd
{

namespace symbolic
{

ExpandVisitor::ExpandVisitor()
{
}

Sum ExpandVisitor::expandedProduct(const Sum& a, const Sum& b)
{
  return a.expandedProduct(b);
}

Sum ExpandVisitor::toExpr(const quotient_map_t& q)
{
  Sum result;
  BOOST_FOREACH(const quotient_map_t::value_type& term, q)
  {
    result += Product::div(term.second, term.first);
  }
  return result;
}

ExpandVisitor::quotient_map_t ExpandVisitor::constructQuotientMap(const Rational& r)
{
  quotient_map_t result;
  result.insert(std::make_pair(Sum(Rational(1)), Sum(r)));
  return result;
}

ExpandVisitor::quotient_map_t ExpandVisitor::reciprocal(const quotient_map_t& q)
{
  // If q has more than one term, this will involve a full expansion.
  Sum nominator;
  Sum denominator(Rational(1));

  BOOST_FOREACH(const quotient_map_t::value_type& qTerm, q)
  {
    nominator = expandedProduct(nominator, qTerm.first) + expandedProduct(denominator, qTerm.second);
    denominator = expandedProduct(denominator, qTerm.first);
  }

  quotient_map_t result;
  // We flip the expanded numerator and denominator
  result.insert(std::make_pair(nominator, denominator));
  return result;
}

void ExpandVisitor::mul(quotient_map_t& q1, const quotient_map_t& q2) const
{
  quotient_map_t result;
  BOOST_FOREACH(const quotient_map_t::value_type& q1Term, q1)
  {
    BOOST_FOREACH(const quotient_map_t::value_type& q2Term, q2)
    {
      const Sum numerator = expandedProduct(q1Term.second, q2Term.second);
      const Sum denominator = expandedProduct(q1Term.first, q2Term.first);
      result[denominator] += numerator;
    }
  }
  std::swap(result, q1);
}

void ExpandVisitor::add(quotient_map_t& q1, const quotient_map_t& q2) const
{
  BOOST_FOREACH(const quotient_map_t::value_type& q2Term, q2)
  {
    q1[q2Term.first] += q2Term.second;
  }
}

void ExpandVisitor::push(const quotient_map_t& q)
{
  stack.push(q);
}

void ExpandVisitor::push(const Sum& s)
{
  const Sum one = Sum(Rational(1));
  quotient_map_t qMap;
  qMap[one] = s;
  push(qMap);
}

void ExpandVisitor::visit(const Sum& s)
{  
  const boost::unordered_map<Expr, quotient_map_t>::const_iterator iter = cache.find(s.clone());

  if (iter == cache.end())
  {
    quotient_map_t qMap = constructQuotientMap(s.getOverall());
    BOOST_FOREACH(const Sum::value_type& term, s)
    {
      term.first.accept(*this);
      quotient_map_t newTerm(stack.top()); stack.pop();
      mul(newTerm, constructQuotientMap(term.second));
      add(qMap, newTerm);
    }

    push(qMap);
    cache.insert(std::make_pair(s.clone(), qMap));
  }
  else
  {
    push(iter->second);
  }
}

void ExpandVisitor::visit(const Product& p)
{
  const boost::unordered_map<Expr, quotient_map_t>::const_iterator iter = cache.find(p.clone());
  if (iter == cache.end())
  {
    quotient_map_t qMap = constructQuotientMap(p.getOverall());
    BOOST_FOREACH(const Product::value_type& term, p)
    {
      term.first.accept(*this);
      quotient_map_t multiplicand(stack.top()); stack.pop();

      if (term.second < 0)
        multiplicand = reciprocal(multiplicand);

      for (int n=0; n < std::abs(term.second); ++n)
        mul(qMap, multiplicand);
    }

    push(qMap);
    cache.insert(std::make_pair(p.clone(), qMap));
  }
  else
  {
    push(iter->second);
  }
}

void ExpandVisitor::visit(const Basic& b)
{
  push(Sum(b));
}

void ExpandVisitor::visit(const Group& g)
{
  push(Sum(g.clone()));
}

Expr ExpandVisitor::getResult() const
{
  assert(stack.size() == 1);
  return toExpr(stack.top()).simplify();
}

}

}
