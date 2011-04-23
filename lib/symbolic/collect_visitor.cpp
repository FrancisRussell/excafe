#include <utility>
#include <cassert>
#include <boost/foreach.hpp>
#include <simple_cfd/symbolic/collect_visitor.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/group.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/basic.hpp>
#include <iostream>

namespace cfd
{

namespace symbolic
{

CollectVisitor::CollectVisitor(const Symbol& s) : symbol(s)
{
  symbolSet.insert(symbol);
}

void CollectVisitor::visit(const Symbol& s)
{
  if (s == symbol)
    stack.push(CollectTerm::poly(s));
  else
    stack.push(CollectTerm::expr(s));
}

void CollectVisitor::visit(const Sum& s)
{
  if (!s.depends(symbolSet))
  {
    stack.push(CollectTerm::expr(s));
  }
  else
  {
    CollectedTerms result = CollectTerm::expr(s.getOverall());

    BOOST_FOREACH(const Sum::value_type& term, s)
    {
      term.first.accept(*this);
      CollectedTerms summand(stack.top()); stack.pop();
      summand *= term.second;
      result += summand;
    }

    stack.push(result);
  }
}

void CollectVisitor::visit(const Product& p)
{
  if (!p.depends(symbolSet))
  {
    stack.push(CollectTerm::expr(p));
  }
  else
  {
    CollectedTerms result = CollectTerm::expr(p.getOverall());

    BOOST_FOREACH(const Product::value_type& term, p)
    {
      if (term.first.depends(symbolSet) && term.second >= 0)
      {
        term.first.accept(*this);
        const CollectedTerms multiplicand(stack.top()); stack.pop();

        for(int i=0; i < term.second; ++i)
          result *= multiplicand;
      }
      else
      {
        result *= CollectTerm::expr(pow(term.first, term.second));
      }
    }

    stack.push(result);
  }
}

void CollectVisitor::visit(const Basic& b)
{
  stack.push(CollectTerm::expr(b));
}

void CollectVisitor::visit(const Group& g)
{
  if (g.depends(symbolSet))
    g.getExpr().accept(*this);
  else
    stack.push(CollectTerm::expr(g));
}

Expr CollectVisitor::getResult() const
{
  assert(stack.size() == 1);
  return stack.top().toExpr();
}

}

}
