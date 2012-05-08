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

namespace cfd
{

namespace symbolic
{

CollectedTerms::CollectedTerms(const Product& product, const Expr& expr)
{
  termMap->insert(std::make_pair(product, expr));
}

CollectedTerms::CollectedTerms()
{
}

CollectedTerms::iterator CollectedTerms::begin()
{
  return termMap->begin();
}

CollectedTerms::iterator CollectedTerms::end()
{
  return termMap->end();
}

CollectedTerms::const_iterator CollectedTerms::begin() const
{
  return termMap->begin();
}

CollectedTerms::const_iterator CollectedTerms::end() const
{
  return termMap->end();
}

CollectedTerms CollectedTerms::poly(const Symbol& s)
{
  return CollectedTerms(Product(s), Rational(1).clone());
}

CollectedTerms CollectedTerms::expr(const Expr& e)
{
  return CollectedTerms(Product::constant(1), e);
}

CollectedTerms& CollectedTerms::operator+=(const CollectedTerms& c)
{
  BOOST_FOREACH(const TermMap::value_type& term, *c.termMap)
    (*termMap)[term.first] += term.second;

  return *this;
}

CollectedTerms& CollectedTerms::operator*=(const CollectedTerms& c)
{
  LazyTermMap resultMap;
  BOOST_FOREACH(const TermMap::value_type& aTerm, *termMap)
  {
    BOOST_FOREACH(const TermMap::value_type& bTerm, *c.termMap)
    {
      const Product product = aTerm.first * bTerm.first;
      (*resultMap)[product] += aTerm.second * bTerm.second;
    }
  }

  std::swap(resultMap, termMap);
  return *this;
}

CollectedTerms& CollectedTerms::operator*=(const Rational& r) 
{
  BOOST_FOREACH(TermMap::value_type& term, *termMap)
  {
    term.second *= r;
  }
  return *this;
}

Expr CollectedTerms::toExpr() const
{
  Sum result;
  BOOST_FOREACH(const TermMap::value_type& term, *termMap)
    result += Product::mul(term.first, term.second);

  return result.simplify();
}


// class CollectVisitor

CollectVisitor::CollectVisitor(const Symbol& s)
{
  symbols.insert(s);
}

CollectVisitor::CollectVisitor(const std::set<Symbol>& _symbols) : symbols(_symbols)
{
}

void CollectVisitor::visit(const Symbol& s)
{
  if (symbols.find(s) != symbols.end())
    stack.push(CollectedTerms::poly(s));
  else
    stack.push(CollectedTerms::expr(s));
}

void CollectVisitor::visit(const Sum& s)
{
  if (!s.depends(symbols))
  {
    stack.push(CollectedTerms::expr(s));
  }
  else
  {
    CollectedTerms result = CollectedTerms::expr(s.getOverall());

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
  if (!p.depends(symbols))
  {
    stack.push(CollectedTerms::expr(p));
  }
  else
  {
    CollectedTerms result = CollectedTerms::expr(p.getOverall());

    BOOST_FOREACH(const Product::value_type& term, p)
    {
      if (term.first.depends(symbols) && term.second >= 0)
      {
        term.first.accept(*this);
        const CollectedTerms multiplicand(stack.top()); stack.pop();

        for(int i=0; i < term.second; ++i)
          result *= multiplicand;
      }
      else
      {
        result *= CollectedTerms::expr(pow(term.first, term.second));
      }
    }

    stack.push(result);
  }
}

void CollectVisitor::visit(const Basic& b)
{
  stack.push(CollectedTerms::expr(b));
}

void CollectVisitor::visit(const Group& g)
{
  if (g.depends(symbols))
    g.getExpr().accept(*this);
  else
    stack.push(CollectedTerms::expr(g));
}

Expr CollectVisitor::getResult() const
{
  assert(stack.size() == 1);
  return stack.top().toExpr();
}

Expr CollectVisitor::getIntegratedResult(const Expr::region_t& region, const unsigned flags) const
{
  assert(stack.size() == 1);
  
  Sum result;
  BOOST_FOREACH(const CollectedTerms::value_type& term, stack.top())
  {
    if (!term.second.depends(symbols))
    {
      result += term.first.integrate(region, flags) * term.second;
    }
    else
    {
      result += Product::mul(term.first, term.second).integrate(region, flags);
    }
  }

  return result.simplify();
}


}

}
