#include <utility>
#include <cassert>
#include <boost/foreach.hpp>
#include <simple_cfd/symbolic/expand_visitor.hpp>
#include <simple_cfd/symbolic/rational.hpp>

namespace cfd
{

namespace symbolic
{

void ExpandVisitor::push(const Expr& e)
{
  const Expr simplified = e.simplify();
  
  if (is_a<Sum>(simplified))
    stack.push(convert_to<Sum>(simplified));
  else
    stack.push(Sum(simplified));
}

void ExpandVisitor::visit(const Sum& s)
{
  Sum reduction(s.getOverall());
  BOOST_FOREACH(const Sum::value_type& term, std::make_pair(s.begin(), s.end()))
  {
    term.first.accept(*this);
    stack.push(Sum(term.second));

    const Sum a = stack.top(); stack.pop();
    const Sum b = stack.top(); stack.pop();
    reduction = reduction + a.expandedProduct(b);
  }
  push(reduction);
}

void ExpandVisitor::visit(const Product& p)
{
  Sum dividend(p.getOverall());
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

  const Product quotient = Product::mul(dividend, Product::pow(divisor, -1));
  push(quotient);
}

void ExpandVisitor::visit(const Basic& b)
{
  push(b);
}

Sum ExpandVisitor::getResult() const
{
  assert(stack.size() == 1);
  return stack.top();
}

}

}
