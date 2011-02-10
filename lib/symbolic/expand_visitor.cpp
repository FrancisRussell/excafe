#include <utility>
#include <cassert>
#include <boost/foreach.hpp>
#include <simple_cfd/symbolic/expand_visitor.hpp>
#include <simple_cfd/symbolic/rational.hpp>

namespace cfd
{

namespace symbolic
{

void ExpandVisitor::visit(const Sum& s)
{
  Sum reduction;
  BOOST_FOREACH(const Sum::value_type& term, std::make_pair(s.begin(), s.end()))
  {
    term.first.accept(*this);
    stack.push(Sum(term.second));

    const Expr a = stack.top(); stack.pop();
    const Expr b = stack.top(); stack.pop();
    reduction = reduction + a*b;
  }

  stack.push(reduction);
}

void ExpandVisitor::visit(const Product& p)
{
  Sum dividend(Rational(1));
  Sum divisor(Rational(1));

  BOOST_FOREACH(const Product::value_type& term, std::make_pair(p.begin(), p.end()))
  {
    term.first.accept(*this);
    const Sum multiplicand(stack.top()); stack.pop();

    Sum& result = (term.second < 0) ? divisor : dividend;
    for(int i=0; i<std::abs(term.second); ++i)
    {
      result = result.expandedProduct(multiplicand);
    }
  }

  const Product quotient(dividend, Product::pow(divisor, -1));
  stack.push(Sum(quotient));
}

void ExpandVisitor::visit(const Basic& b)
{
  stack.push(Sum(b.clone()));
}

Sum ExpandVisitor::getResult() const
{
  assert(stack.size() == 1);
  return stack.top();
}

}

}
