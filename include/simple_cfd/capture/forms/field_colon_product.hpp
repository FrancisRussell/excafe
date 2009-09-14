#ifndef SIMPLE_CFD_CAPTURE_FORMS_COLON_PRODUCT_HPP
#define SIMPLE_CFD_CAPTURE_FORMS_COLON_PRODUCT_HPP

#include <cstddef>
#include "field_expr.hpp"
#include "field_binary_operator.hpp"

namespace cfd
{

namespace detail
{

class FieldColonProduct : public FieldBinaryOperator
{
public:
  FieldColonProduct(FieldExpr::reference_t l, FieldExpr::reference_t r) : FieldBinaryOperator(l, r)
  {
  }

  virtual void accept(FieldVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
