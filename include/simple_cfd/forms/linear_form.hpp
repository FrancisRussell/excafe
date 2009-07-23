#ifndef SIMPLE_CFD_FORMS_LINEAR_FORM_HPP
#define SIMPLE_CFD_FORMS_LINEAR_FORM_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include <simple_cfd/finite_element.hpp>
#include <simple_cfd/fe_vector.hpp>

#include "field.hpp"
#include "basis_field.hpp"
#include "discrete_field.hpp"
#include "tensor_literal.hpp"
#include "facet_normal.hpp"

namespace cfd
{

namespace forms
{

class LinearForm
{
private:
  Field::reference_t field;

public:
  LinearForm(Field* const f) : field(f)
  {
  }

  LinearForm(Field::reference_t f) : field(f)
  {
  }

  template<std::size_t D>
  LinearForm(const FiniteElement<D>& element) : field(new BasisField(element))
  {
  }

  template<std::size_t D>
  LinearForm(const FEVector<D>& element) : field(new DiscreteField(element))
  {
  }

  template<std::size_t D>
  LinearForm(const Tensor<D>& tensor) : field(new TensorLiteral(tensor))
  {
  }

  LinearForm(const facet_normal_tag& tag) : field(new FacetNormal())
  {
  }

  boost::shared_ptr<Field> getField() const
  {
    return field;
  }

  void accept(FieldVisitor& visitor) const
  {
    field->accept(visitor);
  }
};

}

}

#endif
