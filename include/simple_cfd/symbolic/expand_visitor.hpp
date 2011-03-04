#ifndef SIMPLE_CFD_SYMBOLIC_EXPAND_VISITOR_HPP
#define SIMPLE_CFD_SYMBOLIC_EXPAND_VISITOR_HPP

#include <stack>
#include "symbolic_fwd.hpp"
#include "visitor.hpp"
#include "sum.hpp"
#include "product.hpp"
#include "group.hpp"

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
  std::stack<Sum> stack;
  void push(const Expr& e);

public:
  void visit(const Sum& s);
  void visit(const Product& p);
  void visit(const Basic& b);
  void visit(const Group& g);
  Sum getResult() const;
};

}

}

#endif
