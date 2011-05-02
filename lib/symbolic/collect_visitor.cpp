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

// class CollectTerm

CollectTerm& CollectTerm::normalise()
{
  const Rational content = polynomial.findMultiplier();

  if (content == 0)
  {
    polynomial = Sum::constant(0);
    expression = Rational(0);
  }
  else if (content != 1)
  {
    polynomial /= content;
    expression *= content;
  }

  expression = expression.simplify();
  return *this;
}


// class CollectedTerms

void CollectedTerms::addTerm(TermSet& termSet, const CollectTerm& term)
{
  if (term.isZero())
    return;

  // Add for matching polynomial.
  if (addTermTagged<poly_tag>(termSet, term))
    return;

  // Add for matching expression.
  if (addTermTagged<expr_tag>(termSet, term))
    return;

  // If all else fails, create a new term.
  termSet.insert(term);
}

CollectedTerms::CollectedTerms(const CollectTerm& t)
{
  *this += t;
}

CollectedTerms& CollectedTerms::operator+=(const CollectedTerms& c)
{
  BOOST_FOREACH(const TermSet::value_type& term, *c.termSet)
    *this += term;

  return *this;
}

CollectedTerms& CollectedTerms::operator+=(const CollectTerm& t)
{
  addTerm(*termSet, t);
  return *this;
}

CollectedTerms& CollectedTerms::operator*=(const CollectedTerms& c)
{
  LazyTermSet resultSet;
  BOOST_FOREACH(const TermSet::value_type& aTerm, *termSet)
  {
    BOOST_FOREACH(const TermSet::value_type& bTerm, *c.termSet)
    {
      CollectTerm newTerm = aTerm;
      newTerm *= bTerm;
      addTerm(*resultSet, newTerm);
    }
  }

  std::swap(resultSet, termSet);
  return *this;
}

CollectedTerms& CollectedTerms::operator*=(const Rational& r) 
{
  // Since multiplication changes the expression but not the
  // polynomial, we use the polynomial index.

  typedef TermSet::index<poly_tag>::type PolyTermSet;
  PolyTermSet& polyTermSet = termSet->get<poly_tag>();

  for (PolyTermSet::iterator termIter = polyTermSet.begin(); termIter != polyTermSet.end(); ++termIter)
  {
    CollectTerm term = *termIter;
    term *= r;

    const bool replaced = replacePreservingIter(polyTermSet, termIter, term);
    assert(replaced);
  }

  return *this;
}

Expr CollectedTerms::toExpr() const
{
  Sum result;
  BOOST_FOREACH(const TermSet::value_type& aTerm, *termSet)
    result += aTerm.toExpr();

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
    stack.push(CollectTerm::poly(s));
  else
    stack.push(CollectTerm::expr(s));
}

void CollectVisitor::visit(const Sum& s)
{
  if (!s.depends(symbols))
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
  if (!p.depends(symbols))
  {
    stack.push(CollectTerm::expr(p));
  }
  else
  {
    CollectedTerms result = CollectTerm::expr(p.getOverall());

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
  if (g.depends(symbols))
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
