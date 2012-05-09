#ifndef EXCAFE_CAPTURE_FORMS_COLON_PRODUCT_HPP
#define EXCAFE_CAPTURE_FORMS_COLON_PRODUCT_HPP

#include <cstddef>
#include "field_expr.hpp"
#include "field_binary_operator.hpp"

namespace excafe
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
    v.enter(*this);
    getLeft().accept(v);
    getRight().accept(v);
    v.exit(*this);
  }
};

}

}

#endif
