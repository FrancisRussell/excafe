#ifndef SIMPLE_CFD_CAPTURE_FORMS_FIELD_UNARY_OPERATOR_HPP
#define SIMPLE_CFD_CAPTURE_FORMS_FIELD_UNARY_OPERATOR_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "field_expr.hpp"

namespace cfd
{

namespace detail
{

class FieldUnaryOperator : public FieldExpr
{
private:
  FieldExpr::reference_t operand;

public:
  FieldUnaryOperator(FieldExpr::reference_t f) : operand(f)
  {
  }

  FieldExpr::reference_t getOperand() const
  {
    return operand;
  }
};

}

}

#endif
