#ifndef SIMPLE_CFD_FORMS_UNARY_OPERATOR_HPP
#define SIMPLE_CFD_FORMS_UNARY_OPERATOR_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "field.hpp"

namespace cfd
{

namespace forms
{

class UnaryOperator : public Field
{
private:
  Field::reference_t operand;

public:
  UnaryOperator(Field::reference_t f) : operand(f)
  {
  }

  Field::reference_t getOperand() const
  {
    return operand;
  }
};

}

}

#endif
