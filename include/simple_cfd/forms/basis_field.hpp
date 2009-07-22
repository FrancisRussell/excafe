#ifndef SIMPLE_CFD_FORMS_BASIS_FIELD_HPP
#define SIMPLE_CFD_FORMS_BASIS_FIELD_HPP

#include <cstddef>
#include <simple_cfd/finite_element.hpp>
#include "holders.hpp"

namespace cfd
{

namespace forms
{

class BasisField : public Field
{
private:
  FiniteElementHolder element;

public:
  template<std::size_t D>
  BasisField(const FiniteElement<D>& _element) : element(_element)
  {
  }

  std::size_t getRank() const
  {
    return element.getRank();
  }

  virtual void accept(FieldVisitor& v)
  {
    v.visit(*this);
  }

  FiniteElementHolder getElement() const
  {
    return element;
  }
};

}

}
#endif
