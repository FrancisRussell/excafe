#ifndef SIMPLE_CFD_CAPTURE_FORMS_FIELD_SCALAR_HPP
#define SIMPLE_CFD_CAPTURE_FORMS_FIELD_SCALAR_HPP

#include <simple_cfd/capture/fields/scalar.hpp>
#include "field_expr.hpp"

namespace cfd
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
