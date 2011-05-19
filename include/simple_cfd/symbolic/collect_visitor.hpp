#ifndef SIMPLE_CFD_SYMBOLIC_COLLECT_VISITOR_HPP
#define SIMPLE_CFD_SYMBOLIC_COLLECT_VISITOR_HPP

#include <stack>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <boost/unordered_map.hpp>
#include <simple_cfd/exception.hpp>
#include "symbolic_fwd.hpp"
#include "visitor.hpp"
#include "sum.hpp"
#include "expr.hpp"
#include "product.hpp"
#include "group.hpp"
#include "symbol.hpp"
#include <simple_cfd/util/lazy_copy.hpp>

namespace cfd
{

namespace symbolic
{

class CollectedTerms
{
private:
  typedef boost::unordered_map<Sum, Expr> TermMap;
  typedef util::LazyCopy<TermMap> LazyTermMap;

  LazyTermMap termMap;

  CollectedTerms(const Sum& sum, const Expr& expr);

public:
  static CollectedTerms poly(const Symbol& s);
  static CollectedTerms expr(const Expr& e);

  CollectedTerms();
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
};

}

}

#endif
