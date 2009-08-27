#ifndef SIMPLE_CFD_FORMS_DIVERGENCE_HPP
#define SIMPLE_CFD_FORMS_DIVERGENCE_HPP

#include <cstddef>
#include <cassert>
#include "field_expr.hpp"
#include "field_unary_operator.hpp"

namespace cfd
{

namespace detail
{

class FieldDivergence : public FieldUnaryOperator
{
public:
  FieldDivergence(FieldExpr::reference_t f) : FieldUnaryOperator(f)
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
