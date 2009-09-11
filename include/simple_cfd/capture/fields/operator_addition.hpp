#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ADDITION_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ADDITION_HPP

#include <cassert>
#include "operator_expr.hpp"
#include "discrete_expr_visitor.hpp"

namespace cfd
{

namespace detail
{

class OperatorAddition : public OperatorExpr
{
private:
  OperatorExpr::expr_ptr left;
  OperatorExpr::expr_ptr right;

public:
  OperatorAddition(const OperatorExpr::expr_ptr& _left, const OperatorExpr::expr_ptr& _right) :
    left(_left), right(_right)
  {
  }

  virtual void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getTrialSpace() const
  {
    assert(left->getTrialSpace() == right->getTrialSpace());
    return left->getTrialSpace();
  }

  virtual FunctionSpaceExpr::expr_ptr getTestSpace() const
  {
    assert(left->getTestSpace() == right->getTestSpace());
    return left->getTestSpace();
  }
};

}

}

#endif
