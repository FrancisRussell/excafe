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
  BOOST_FOREACH(const Sum::value_type& term, std::make_pair(s.begin(), s.end()))
  {
    term.first.accept(*this);
    stack.push(Sum::multiplier(term.second));

    const Sum a = stack.top(); stack.pop();
    const Sum b = stack.top(); stack.pop();
    stack.push(a*b);
  }

  Sum reduction;
  assert(stack.size() >= s.nops());
  for(std::size_t i=0; i<s.nops(); ++i)
  {
    reduction = reduction + stack.top();
    stack.pop();
  }

  stack.push(reduction);
}

void ExpandVisitor::accept(const Product& p)
{
  BOOST_FOREACH(const Product::value_type& term, std::make_pair(p.begin(), p.end()))
  {
    term.first.accept(*this);

    const Sum multiplicand(stack.top()); stack.pop();
    Sum result = Sum::multiplier(1);

    for(int i=0; i<term.second; ++i)
      result = result.expandedProduct(multiplicand);

    stack.push(result);
  }

  Sum reduction = Sum::multiplier(1);
  assert(stack.size() >= p.nops());
  for(std::size_t i=0; i<p.nops(); ++i)
  {
    reduction = reduction.expandedProduct(stack.top());
    stack.pop();
  }

  stack.push(reduction);
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
