#ifndef SIMPLE_CFD_FORMS_GRADIENT_HPP
#define SIMPLE_CFD_FORMS_GRADIENT_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "unary_operator.hpp"
#include "field.hpp"

namespace cfd
{

namespace forms
{

class Gradient : public UnaryOperator
{
public:
  Gradient(Field::reference_t f) : UnaryOperator(f)
  {
  }

  std::size_t getRank() const
  {
    return getOperand()->getRank() + 1;
  }

  std::size_t getDimension() const
  {
    return getOperand()->getDimension();
  }

  virtual void accept(FieldVisitor& v)
  {
    v.enter(*this);
    getOperand()->accept(v);
    v.exit(*this);
  }
};

}

}

#endif
