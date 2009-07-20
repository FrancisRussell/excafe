#ifndef SIMPLE_CFD_FORMS_BINARY_OPERATOR_HPP
#define SIMPLE_CFD_FORMS_BINARY_OPERATOR_HPP

#include <cstddef>
#include "field.hpp"

namespace cfd
{

namespace forms
{

class BinaryOperator : public Field
{
private:
  Field::reference_t left;
  Field::reference_t right;

public:
  BinaryOperator(Field::reference_t l, Field::reference_t r) : left(l), right(r)
  {
  }

  Field::reference_t getLeft() const
  {
    return left;
  }

  Field::reference_t getRight() const
  {
    return right;
  }
};

}

}

#endif
