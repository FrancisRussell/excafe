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

namespace cfd
{

namespace forms
{

class LinearForm
{
private:
  boost::shared_ptr<Field> field;

public:
  LinearForm(Field* const f) : field(f)
  {
  }

  LinearForm(boost::shared_ptr<Field> f) : field(f)
  {
  }

  template<std::size_t D>
  LinearForm(const FiniteElement<D>& element) : field(new BasisField(element))
  {
  }

  template<typename C>
  LinearForm(const FEVector<C>& element) : field(new DiscreteField(element))
  {
  }

  template<std::size_t D>
  LinearForm(const Tensor<D>& tensor) : field(new TensorLiteral(tensor))
  {
  }

  boost::shared_ptr<Field> getField() const
  {
    return field;
  }
};

}

}

#endif
