#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ADDITION_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ADDITION_HPP

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
  OperatorAddition(const OperatorExpr::expr_ptr _left, const OperatorExpr::expr_ptr _right) :
    left(_left), right(_right)
  {
  }

  virtual void accept(DiscreteExprVisitor& v)
  {
    v.enter(*this);
    left->accept(v);
    right->accept(v);
    v.exit(*this);
  }
};

}

}

#endif
