#ifndef SIMPLE_CFD_CAPTURE_FORMS_FIELD_INNER_PRODUCT_HPP
#define SIMPLE_CFD_CAPTURE_FORMS_FIELD_INNER_PRODUCT_HPP

#include <cstddef>
#include "field_expr.hpp"
#include "field_binary_operator.hpp"

namespace cfd
{

namespace detail
{

class FieldInnerProduct : public FieldBinaryOperator
{
public:
  FieldInnerProduct(FieldExpr::reference_t l, FieldExpr::reference_t r) : FieldBinaryOperator(l, r)
  {
  }

  virtual void accept(FieldVisitor& v)
  {
    v.enter(*this);
    getLeft()->accept(v);
    getRight()->accept(v);
    v.exit(*this);
  }
};

}

}
#endif
