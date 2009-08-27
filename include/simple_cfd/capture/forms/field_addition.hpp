#ifndef SIMPLE_CFD_CAPTURE_FORMS_FIELD_ADDITION_HPP
#define SIMPLE_CFD_CAPTURE_FORMS_FIELD_ADDITION_HPP

#include <cstddef>
#include <cassert>
#include "field_expr.hpp"
#include "field_binary_operator.hpp"

namespace cfd
{

namespace detail
{

class FieldAddition : public FieldBinaryOperator
{
public:
  FieldAddition(FieldExpr::reference_t l, FieldExpr::reference_t r) : FieldBinaryOperator(l, r)
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
