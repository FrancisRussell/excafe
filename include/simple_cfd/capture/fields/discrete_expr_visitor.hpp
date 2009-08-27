#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_VISITOR_HPP

#include "fields_fwd.hpp"

namespace cfd
{

namespace detail
{

class DiscreteExprVisitor
{
public:
  // Discrete field related
  virtual void visit(DiscreteFieldUndefined& u) = 0;
  virtual void visit(DiscreteFieldZero& z) = 0;
  virtual void visit(DiscreteFieldPersistent& p) = 0;

  // Discrete operator related
  virtual void enter(OperatorApplication& a) = 0;
  virtual void exit(OperatorApplication& a) = 0;
  virtual void visit(OperatorAssembly& a) = 0;
  virtual void visit(OperatorUndefined& u) = 0;

  // Scalar related
  virtual void enter(ScalarBinaryOperator& o) = 0;
  virtual void exit(ScalarBinaryOperator& o) = 0;
  virtual void visit(ScalarLiteral& l) = 0;

  virtual ~DiscreteExprVisitor() {}
};

}

}


#endif
