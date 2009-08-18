#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELD_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELD_HPP

#include "field_expr.hpp"
#include "field_empty.hpp"

namespace cfd
{

class Field
{
private:
  typedef detail::FieldExpr::expr_ptr expr_ptr;
  expr_ptr expr;

public:
  Field() : expr(new detail::FieldEmpty())
  {
  }

  Field(detail::FieldExpr* const _expr) : expr(_expr)
  {
  }
};

}

#endif
