#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_VISITOR_HPP

#include "fields_fwd.hpp"
#include "discrete_traits.hpp"

namespace cfd
{

namespace detail
{

class DiscreteExprVisitor
{
public:
  // Discrete field related
  virtual void enter(DiscreteFieldElementWise& p) = 0;
  virtual void exit(DiscreteFieldElementWise& p) = 0;
  virtual void enter(DiscreteFieldTwoNorm& p) = 0;
  virtual void exit(DiscreteFieldTwoNorm& p) = 0;
  virtual void visit(DiscreteFieldUndefined& u) = 0;
  virtual void visit(DiscreteFieldZero& z) = 0;
  virtual void visit(DiscreteFieldPersistent& p) = 0;

  // Discrete operator related
  virtual void enter(OperatorApplication& a) = 0;
  virtual void exit(OperatorApplication& a) = 0;
  virtual void enter(OperatorAddition& u) = 0;
  virtual void exit(OperatorAddition& u) = 0;
  virtual void visit(OperatorAssembly& a) = 0;
  virtual void visit(OperatorUndefined& u) = 0;

  // Scalar related
  virtual void enter(ScalarBinaryOperator& o) = 0;
  virtual void exit(ScalarBinaryOperator& o) = 0;
  virtual void visit(ScalarLiteral& l) = 0;

  // Temporal related
  virtual void visit(DiscreteObjectIndexed<discrete_scalar_tag>& s) = 0;
  virtual void visit(DiscreteObjectIndexed<discrete_field_tag>& s) = 0;
  virtual void visit(DiscreteObjectIndexed<discrete_operator_tag>& s) = 0;

  // Solve related
  virtual void enter(LinearSolve& s) = 0;
  virtual void exit(LinearSolve& s) = 0;

  virtual ~DiscreteExprVisitor() {}
};

}

}


#endif
