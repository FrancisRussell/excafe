#ifndef EXCAFE_SYMBOLIC_EXPAND_VISITOR_HPP
#define EXCAFE_SYMBOLIC_EXPAND_VISITOR_HPP

#include <stack>
#include <boost/unordered_map.hpp>
#include "symbolic_fwd.hpp"
#include "visitor.hpp"
#include "sum.hpp"
#include "product.hpp"
#include "group.hpp"
#include "symbol.hpp"

namespace excafe
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

  static Sum toExpr(const quotient_map_t& q, bool distribute);
  static quotient_map_t constructQuotientMap(const Rational& r);
  static quotient_map_t reciprocal(const quotient_map_t& q);
  static Sum expandedProduct(const Sum& a, const Sum& b);

  std::stack<quotient_map_t> stack;

  void add(quotient_map_t& q1, const quotient_map_t& q2) const;
  void mul(quotient_map_t& q1, const quotient_map_t& q2) const;
  void mul(quotient_map_t& q1, const Rational& r) const;
  void push(const quotient_map_t& e);
  void push(const Sum& s);

public:
  ExpandVisitor();
  void visit(const Sum& s);
  void visit(const Product& p);
  void visit(const Basic& b);
  void visit(const Group& g);
  Expr getResult(bool distribute) const;
};

}

}

#endif
