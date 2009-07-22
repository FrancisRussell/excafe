#ifndef SIMPLE_CFD_FORMS_DIVERGENCE_HPP
#define SIMPLE_CFD_FORMS_DIVERGENCE_HPP

#include <cstddef>
#include <cassert>
#include "unary_operator.hpp"

namespace cfd
{

namespace forms
{

class Divergence : public UnaryOperator
{
public:
  Divergence(Field::reference_t f) : UnaryOperator(f)
  {
    assert(getOperand()->getRank() > 0);
  }

  std::size_t getRank() const
  {
    return getOperand()->getRank() - 1;
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
