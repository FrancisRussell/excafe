#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_UNDEFINED_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_UNDEFINED_HPP

#include "discrete_expr_visitor.hpp"
#include "operator_expr.hpp"

namespace cfd
{

namespace detail
{

class OperatorUndefined : public OperatorExpr
{
public:
  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
