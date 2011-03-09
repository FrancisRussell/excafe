#ifndef SIMPLE_CFD_SYMBOLIC_EXPAND_VISITOR_HPP
#define SIMPLE_CFD_SYMBOLIC_EXPAND_VISITOR_HPP

#include <stack>
#include <boost/optional.hpp>
#include <boost/tr1/unordered_map.hpp>
#include "symbolic_fwd.hpp"
#include "visitor.hpp"
#include "sum.hpp"
#include "product.hpp"
#include "group.hpp"
#include "symbol.hpp"

namespace cfd
{

namespace symbolic
{

class ExpandVisitor : public Visitor, 
                      public Sum::Visitor,
                      public Product::Visitor,
                      public Group::Visitor
{
private:
  typedef boost::unordered_map<Sum, Sum> quotient_map_t;
  static Sum toExpr(const quotient_map_t& q);
  static quotient_map_t constructQuotientMap(const Rational& r);
  quotient_map_t reciprocal(const quotient_map_t& q) const;
  void add(quotient_map_t& q1, const quotient_map_t& q2) const;
  void mul(quotient_map_t& q1, const quotient_map_t& q2) const;

  boost::optional<Symbol> symbol;
  std::stack<quotient_map_t> stack;
  void push(const quotient_map_t& e);
  void push(const Sum& s);

public:
  ExpandVisitor();
  ExpandVisitor(const Symbol& s);

  void visit(const Sum& s);
  void visit(const Product& p);
  void visit(const Basic& b);
  void visit(const Group& g);
  Expr getResult() const;
};

}

}

#endif
