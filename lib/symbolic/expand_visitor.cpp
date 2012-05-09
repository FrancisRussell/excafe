#include <utility>
#include <cassert>
#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <excafe/symbolic/expand_visitor.hpp>
#include <excafe/symbolic/rational.hpp>
#include <excafe/symbolic/group.hpp>
#include <excafe/symbolic/product.hpp>
#include <excafe/symbolic/sum.hpp>
#include <excafe/symbolic/expr.hpp>
#include <excafe/symbolic/basic.hpp>

namespace excafe
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

Sum ExpandVisitor::toExpr(const quotient_map_t& q, const bool distribute)
{
  Sum result;
  BOOST_FOREACH(const quotient_map_t::value_type& term, q)
  {
    if (distribute)
      result += expandedProduct(term.second, Sum(Product::div(Rational(1), term.first)));
    else
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

void ExpandVisitor::mul(quotient_map_t& q1, const Rational& r) const
{
  BOOST_FOREACH(quotient_map_t::value_type& q1Term, q1)
    q1Term.second *= r;
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
  quotient_map_t qMap = constructQuotientMap(s.getOverall());
  BOOST_FOREACH(const Sum::value_type& term, s)
  {
    term.first.accept(*this);
    quotient_map_t newTerm(stack.top()); stack.pop();
    mul(newTerm, term.second);
    add(qMap, newTerm);
  }

  push(qMap);
}

void ExpandVisitor::visit(const Product& p)
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
}

void ExpandVisitor::visit(const Basic& b)
{
  push(Sum(b));
}

void ExpandVisitor::visit(const Group& g)
{
  push(Sum(g.clone()));
}

Expr ExpandVisitor::getResult(const bool distribute) const
{
  assert(stack.size() == 1);
  return toExpr(stack.top(), distribute).simplify();
}

}

}
