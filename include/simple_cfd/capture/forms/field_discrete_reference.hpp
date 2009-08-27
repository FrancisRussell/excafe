#ifndef SIMPLE_CFD_CAPTURE_FORMS_FIELD_DISCRETE_REFERENCE_HPP
#define SIMPLE_CFD_CAPTURE_FORMS_FIELD_DISCRETE_REFERENCE_HPP

#include <cstddef>
#include "field_expr.hpp"
#include <simple_cfd/capture/fields/field.hpp>

namespace cfd
{

namespace detail
{

class FieldDiscreteReference : public FieldExpr
{
private:
  Field discreteField;

public:
  FieldDiscreteReference(const Field& field) : discreteField(field)
  {
  }

  virtual void accept(FieldVisitor& v)
  {
    v.visit(*this);
  }

  Field getDiscreteField() const
  {
    return discreteField;
  }
};

}

}

#endif
