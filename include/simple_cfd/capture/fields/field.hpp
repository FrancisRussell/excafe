#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELD_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELD_HPP

#include "fields_fwd.hpp"
#include "discrete_field_expr.hpp"

namespace cfd
{

class Field
{
public:
  typedef detail::DiscreteFieldExpr::expr_ptr expr_ptr;

private:
  expr_ptr expr;

public:
  Field();
  Field(const FunctionSpace& functionSpace);
  Field(detail::DiscreteFieldExpr* const _expr);
  expr_ptr getExpr() const;
};

}

#endif
