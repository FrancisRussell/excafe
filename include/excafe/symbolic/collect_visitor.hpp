#ifndef EXCAFE_SYMBOLIC_COLLECT_VISITOR_HPP
#define EXCAFE_SYMBOLIC_COLLECT_VISITOR_HPP

#include <stack>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <boost/unordered_map.hpp>
#include <excafe/exception.hpp>
#include "symbolic_fwd.hpp"
#include "visitor.hpp"
#include "expr.hpp"
#include "product.hpp"
#include "sum.hpp"
#include "group.hpp"
#include "symbol.hpp"
#include <excafe/util/lazy_copy.hpp>

namespace excafe
{

namespace symbolic
{

class CollectedTerms
{
private:
  typedef boost::unordered_map<Product, Expr> TermMap;
  typedef util::LazyCopy<TermMap> LazyTermMap;

  LazyTermMap termMap;

  CollectedTerms(const Product& sum, const Expr& expr);

public:
  typedef TermMap::value_type value_type;
  typedef TermMap::iterator iterator;
  typedef TermMap::const_iterator const_iterator;

  static CollectedTerms poly(const Symbol& s);
  static CollectedTerms expr(const Expr& e);

  CollectedTerms();
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;
  CollectedTerms& operator+=(const CollectedTerms& c);
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
  std::set<Symbol> symbols;
  std::stack<CollectedTerms> stack;

public:
  CollectVisitor(const Symbol& s);
  CollectVisitor(const std::set<Symbol>& s);

  void visit(const Basic& b);
  void visit(const Group& g);
  void visit(const Symbol& s);
  void visit(const Sum& s);
  void visit(const Product& p);
  Expr getResult() const;
  Expr getIntegratedResult(const Expr::region_t& region, const unsigned flags) const;
};

}

}

#endif
