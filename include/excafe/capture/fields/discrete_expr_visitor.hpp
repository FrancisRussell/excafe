#ifndef EXCAFE_CAPTURE_FIELDS_DISCRETE_EXPR_VISITOR_HPP
#define EXCAFE_CAPTURE_FIELDS_DISCRETE_EXPR_VISITOR_HPP

#include "fields_fwd.hpp"

namespace excafe
{

namespace detail
{

class DiscreteExprVisitor
{
public:
  // Discrete field related
  virtual void visit(DiscreteFieldElementWise& p) = 0;
  virtual void visit(DiscreteFieldTwoNorm& p) = 0;
  virtual void visit(DiscreteFieldProjection& p) = 0;
  virtual void visit(DiscreteFieldUndefined& u) = 0;
  virtual void visit(DiscreteFieldZero& z) = 0;
  virtual void visit(DiscreteFieldPersistent& p) = 0;
  virtual void visit(DiscreteFieldApplyBC& a) = 0;

  // Discrete operator related
  virtual void visit(OperatorApplication& a) = 0;
  virtual void visit(OperatorAddition& u) = 0;
  virtual void visit(OperatorAssembly& a) = 0;
  virtual void visit(OperatorUndefined& u) = 0;
  virtual void visit(OperatorApplyBC& a) = 0;

  // Scalar related
  virtual void visit(ScalarBinaryOperator& o) = 0;
  virtual void visit(ScalarLiteral& l) = 0;
  virtual void visit(ScalarUndefined& l) = 0;

  // Temporal related
  virtual void visit(DiscreteIndexedScalar& s) = 0;
  virtual void visit(DiscreteIndexedField& s) = 0;
  virtual void visit(DiscreteIndexedOperator& s) = 0;

  // Solve related
  virtual void visit(LinearSolve& s) = 0;

  virtual ~DiscreteExprVisitor() {}
};

}

}


#endif
