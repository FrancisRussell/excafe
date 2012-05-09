#ifndef EXCAFE_FORMS_CAPTURE_FIELD_BINARY_OPERATOR_HPP
#define EXCAFE_FORMS_CAPTURE_FIELD_BINARY_OPERATOR_HPP

#include <cstddef>
#include "field_expr.hpp"

namespace excafe
{

namespace detail
{

class FieldBinaryOperator : public FieldExpr
{
private:
  FieldExpr::reference_t left;
  FieldExpr::reference_t right;

public:
  FieldBinaryOperator(FieldExpr::reference_t l, FieldExpr::reference_t r) : left(l), right(r)
  {
  }

  FieldExpr& getLeft() const
  {
    return *left;
  }

  FieldExpr& getRight() const
  {
    return *right;
  }
};

}

}

#endif
