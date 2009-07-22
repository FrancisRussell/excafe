#ifndef SIMPLE_CFD_FORMS_DISCRETE_FIELD_HPP
#define SIMPLE_CFD_FORMS_DISCRETE_FIELD_HPP

#include <cstddef>
#include "field.hpp"
#include "holders.hpp"
#include <simple_cfd/fe_vector.hpp>

namespace cfd
{

namespace forms
{

class DiscreteField : public Field
{
private:
 FEVectorHolder vector;

public:
  template<std::size_t D>
  DiscreteField(const FEVector<D>& v) : vector(v)
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

  FEVectorHolder getVector() const
  {
    return vector;
  }
};

}

}

#endif
