#ifndef SIMPLE_CFD_SYMBOLIC_EXPAND_VISITOR_HPP
#define SIMPLE_CFD_SYMBOLIC_EXPAND_VISITOR_HPP

#include <stack>
#include "symbolic_fwd.hpp"
#include "visitor.hpp"
#include "sum.hpp"
#include "product.hpp"

namespace cfd
{

namespace symbolic
{

class ExpandVisitor : public Visitor, 
                      public Sum::Visitor,
                      public Product::Visitor
{
private:
  std::stack<Sum> stack;

public:
  void accept(const Sum& s);
  void accept(const Product& p);
  void accept(const Basic& b);
  Sum getResult() const;
};

}

}

#endif
