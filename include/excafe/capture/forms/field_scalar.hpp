#ifndef EXCAFE_CAPTURE_FORMS_FIELD_SCALAR_HPP
#define EXCAFE_CAPTURE_FORMS_FIELD_SCALAR_HPP

#include <excafe/capture/fields/scalar.hpp>
#include "field_expr.hpp"

namespace excafe
{

namespace detail
{

class FieldScalar : public FieldExpr
{
private:
  Scalar scalar;

public:
  FieldScalar(const Scalar& s) : scalar(s)
  {
  }

  virtual void accept(FieldVisitor& v)
  {
    v.visit(*this);
  }

  Scalar getValue() const
  {
    return scalar;
  }
};

}

}

#endif
