#ifndef SIMPLE_CFD_CAPTURE_FIELDS_LINEAR_SOLVE_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_LINEAR_SOLVE_HPP

#include <cassert>
#include "discrete_expr_visitor.hpp"
#include "discrete_field_expr.hpp"

namespace cfd
{

namespace detail
{

class LinearSolve : public DiscreteFieldExpr
{
private:
  OperatorExpr::expr_ptr operation;
  DiscreteFieldExpr::expr_ptr operand;

public:
  LinearSolve(const OperatorExpr::expr_ptr& _operation, const DiscreteFieldExpr::expr_ptr& _operand) :
    operation(_operation), operand(_operand)
  {
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    assert(operation->getTestSpace() == operand->getFunctionSpace());
    return operation->getTrialSpace();
  }
};

}

}

#endif
