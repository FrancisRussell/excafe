#ifndef SIMPLE_CFD_FORMS_DISCRETE_FIELD_HPP
#define SIMPLE_CFD_FORMS_DISCRETE_FIELD_HPP

#include <cstddef>
#include "field.hpp"
#include "holders.hpp"
#include <simple_cfd/discrete_field.hpp>

namespace cfd
{

namespace forms
{

class DiscreteFieldReference : public Field
{
private:
 DiscreteFieldHolder vector;

public:
  template<std::size_t D>
  DiscreteFieldReference(const DiscreteField<D>& v) : vector(v)
  {
  }

  std::size_t getRank() const
  {
    return vector.getRank();
  }

  virtual void accept(FieldVisitor& v)
  {
    v.visit(*this);
  }

  DiscreteFieldHolder getVector() const
  {
    return vector;
  }
};

}

}

#endif
