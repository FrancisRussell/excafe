#ifndef SIMPLE_CFD_SYMBOLIC_COLLECT_VISITOR_HPP
#define SIMPLE_CFD_SYMBOLIC_COLLECT_VISITOR_HPP

#include <stack>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <boost/ref.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/indexed_by.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/tag.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <simple_cfd/exception.hpp>
#include "symbolic_fwd.hpp"
#include "visitor.hpp"
#include "sum.hpp"
#include "expr.hpp"
#include "product.hpp"
#include "group.hpp"
#include "symbol.hpp"

namespace cfd
{

namespace symbolic
{

struct poly_tag {};
struct expr_tag {};

class CollectTerm
{
private:
  Sum  polynomial;
  Expr expression;

  CollectTerm(const Sum& _polynomial, const Expr& _expression) :
    polynomial(_polynomial), expression(_expression)
  {
    normalise();
  }

  CollectTerm& normalise();

public:
  static CollectTerm poly(const Symbol& s) 
  {
    return CollectTerm(Sum(s), Rational(1));
  }

  static CollectTerm expr(const Expr& e)
  {
    return CollectTerm(Sum::constant(1), e);
  }

  const Sum& poly() const
  {
    return polynomial;
  }

  const Expr& expr() const
  {
    return expression;
  }

  bool isZero() const
  {
    return expression == Rational(0);
  }

  CollectTerm& operator+=(const CollectTerm& c)
  {
    if (polynomial == c.polynomial)
    {
      expression += c.expression;
    }
    else if (expression == c.expression)
    {
      polynomial += c.polynomial;
    }
    else
    {
      CFD_EXCEPTION("Can only add CollectTerms when one of the products match.");
    }

    return normalise();
  }

  CollectTerm& operator*=(const CollectTerm& c)
  {
    polynomial = polynomial.expandedProduct(c.polynomial);
    expression *= c.expression;

    return normalise();
  }

  CollectTerm& operator*=(const Rational& r)
  {
    if (r != 1)
      expression *= r;

    return normalise();
  }

  Expr toExpr() const
  {
    return Product::mul(polynomial, expression).simplify();
  }
};

class CollectedTerms
{
private:
  typedef boost::multi_index_container<
    CollectTerm,
    boost::multi_index::indexed_by<
      boost::multi_index::hashed_non_unique<
        boost::multi_index::tag<poly_tag>, 
        BOOST_MULTI_INDEX_CONST_MEM_FUN(CollectTerm, const Sum&, poly)>,
      boost::multi_index::hashed_non_unique<
        boost::multi_index::tag<expr_tag>,
        BOOST_MULTI_INDEX_CONST_MEM_FUN(CollectTerm, const Expr&, expr)>
    >
  > TermSet;

  typedef util::LazyCopy<TermSet> LazyTermSet;

  LazyTermSet termSet;

  template<typename Tag>
  static bool addTermTagged(TermSet& termSet, const CollectTerm& term)
  {
    typedef typename TermSet::index<Tag>::type TaggedTermSet;
    typename TaggedTermSet::key_from_value keyExtractor;
    TaggedTermSet& taggedTermSet = termSet.get<Tag>();

    const typename TaggedTermSet::iterator termIter = taggedTermSet.find(keyExtractor(term));
    if (termIter != taggedTermSet.end())
    {
      CollectTerm newTerm = *termIter;
      newTerm += term;

      if (newTerm.isZero())
      {
        taggedTermSet.erase(termIter);
      }
      else
      {
        const bool replaced = taggedTermSet.replace(termIter, newTerm);
        assert(replaced);
      }

      return true;
    }
    else
    {
      return false;
    }
  }

  template<typename Set>
  static bool replacePreservingIter(Set& set, 
                                    const typename Set::iterator& iterator, 
                                    const typename Set::value_type& value)
  {
    // Asserts that when we replace a value, we replace it with one with the same key. Otherwise, we might
    // encounter the same key again, if using a sorted or hashed data structure.
    typename Set::key_from_value keyExtractor;
    assert(keyExtractor(*iterator) == keyExtractor(value));

    return set.replace(iterator, value);
  }

  static void addTerm(TermSet& termSet, const CollectTerm& term);

public:
  CollectedTerms()
  {
  }

  CollectedTerms(const CollectTerm& t);
  CollectedTerms& operator+=(const CollectedTerms& c);
  CollectedTerms& operator+=(const CollectTerm& t);
  CollectedTerms& operator*=(const CollectedTerms& c);
  CollectedTerms& operator*=(const Rational& r);
  Expr toExpr() const;
};

class CollectVisitor : public Visitor, 
                       public Symbol::Visitor,
                       public Sum::Visitor,
                       public Product::Visitor,
                       public Group::Visitor
{
private:
  Symbol symbol;
  std::set<Symbol> symbolSet;
  std::stack<CollectedTerms> stack;

public:
  CollectVisitor(const Symbol& s);

  void visit(const Basic& b);
  void visit(const Group& g);
  void visit(const Symbol& s);
  void visit(const Sum& s);
  void visit(const Product& p);
  Expr getResult() const;
};

}

}

#endif
