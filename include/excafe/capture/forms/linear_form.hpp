#ifndef EXCAFE_CAPTURE_FORMS_LINEAR_FORM_HPP
#define EXCAFE_CAPTURE_FORMS_LINEAR_FORM_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>

#include <excafe/capture/fields/element.hpp>
#include <excafe/capture/fields/field.hpp>
#include <excafe/capture/fields/named_field.hpp>
#include <excafe/capture/fields/scalar.hpp>
#include <excafe/capture/fields/indexed_value_helper.hpp>
#include <excafe/capture/fields/discrete_traits.hpp>

#include "field_expr.hpp"
#include "field_basis.hpp"
#include "field_discrete_reference.hpp"
#include "field_scalar.hpp"
#include "field_visitor.hpp"
#include "facet_normal.hpp"

namespace excafe
{

namespace forms
{

class LinearForm
{
private:
  detail::FieldExpr::reference_t field;

public:
  class facet_normal_tag {};

  LinearForm(detail::FieldExpr* const f) : field(f)
  {
  }

  LinearForm(const detail::FieldExpr::reference_t f) : field(f)
  {
  }

  LinearForm(const Element& element) : field(new detail::FieldBasis(element))
  {
  }

  LinearForm(const Field& f) : field(new detail::FieldDiscreteReference(f))
  {
  }

  LinearForm(const detail::IndexedValueHelper<detail::discrete_field_tag>& f) : 
    field(new detail::FieldDiscreteReference(f))
  {
  }

  LinearForm(const NamedField& f) : field(new detail::FieldDiscreteReference(f))
  {
  }

  LinearForm(const Scalar& s) : field(new detail::FieldScalar(s))
  {
  }

  LinearForm(const detail::IndexedValueHelper<detail::discrete_scalar_tag>& s) : 
    field(new detail::FieldScalar(s))
  {
  }

  LinearForm(const facet_normal_tag& tag) : field(new detail::FacetNormal())
  {
  }

  detail::FieldExpr::reference_t getField() const
  {
    return field;
  }

  void accept(detail::FieldVisitor& visitor) const
  {
    field->accept(visitor);
  }
};

}

}

#endif
