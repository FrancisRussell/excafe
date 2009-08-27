#ifndef SIMPLE_CFD_CAPTURE_FORMS_GRADIENT_HPP
#define SIMPLE_CFD_CAPTURE_FORMS_GRADIENT_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "field_unary_operator.hpp"
#include "field_expr.hpp"

namespace cfd
{

namespace detail
{

class FieldGradient : public FieldUnaryOperator
{
public:
  FieldGradient(FieldExpr::reference_t f) : FieldUnaryOperator(f)
  {
  }

  virtual void accept(FieldVisitor& v)
  {
    v.enter(*this);
    getOperand()->accept(v);
    v.exit(*this);
  }
};

}

}

#endif
