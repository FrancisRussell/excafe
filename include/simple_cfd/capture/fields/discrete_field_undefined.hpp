#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_EMPTY_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_EMPTY_HPP

#include "discrete_expr_visitor.hpp"
#include "discrete_field_expr.hpp"
#include "function_space_undefined.hpp"

namespace cfd
{

namespace detail
{

class DiscreteFieldUndefined : public DiscreteFieldExpr
{
public:
  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    return FunctionSpaceExpr::expr_ptr(new FunctionSpaceUndefined());
  }
};

}

}

#endif
