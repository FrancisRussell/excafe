#ifndef EXCAFE_CAPTURE_FORMS_FIELD_OUTER_PRODUCT_HPP
#define EXCAFE_CAPTURE_FORMS_FIELD_OUTER_PRODUCT_HPP

#include <cstddef>
#include "field_expr.hpp"
#include "field_binary_operator.hpp"

namespace excafe
{

namespace detail
{

class FieldOuterProduct : public FieldBinaryOperator
{
public:
  FieldOuterProduct(FieldExpr::reference_t l, FieldExpr::reference_t r) : FieldBinaryOperator(l, r)
  {
  }

  virtual void accept(FieldVisitor& v)
  {
    v.enter(*this);
    getLeft().accept(v);
    getRight().accept(v);
    v.exit(*this);
  }
};

}

}

#endif
