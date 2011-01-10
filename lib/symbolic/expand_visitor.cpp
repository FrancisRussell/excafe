#include <utility>
#include <cassert>
#include <boost/foreach.hpp>
#include <simple_cfd/symbolic/expand_visitor.hpp>

namespace cfd
{

namespace symbolic
{

void ExpandVisitor::accept(const Sum& s)
{
  Sum reduction;
  BOOST_FOREACH(const Sum::value_type& term, std::make_pair(s.begin(), s.end()))
  {
    term.first.accept(*this);
    stack.push(Sum::multiplier(term.second));

    const Sum a = stack.top(); stack.pop();
    const Sum b = stack.top(); stack.pop();
    reduction = reduction + a*b;
  }

  stack.push(reduction);
}

void ExpandVisitor::accept(const Product& p)
{
  Sum dividend = Sum::multiplier(1);
  Sum divisor = Sum::multiplier(1);

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

void ExpandVisitor::accept(const Basic& b)
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
