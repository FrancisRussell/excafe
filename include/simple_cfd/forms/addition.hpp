#ifndef SIMPLE_CFD_FORMS_ADDITION_HPP
#define SIMPLE_CFD_FORMS_ADDITION_HPP

#include <cstddef>
#include <cassert>
#include "binary_operator.hpp"

namespace cfd
{

namespace forms
{

class Addition : public BinaryOperator
{
public:
  Addition(Field::reference_t l, Field::reference_t r) : BinaryOperator(l, r)
  {
    assert(l->getRank() == r->getRank());
    assert(l->getDimension() == r->getDimension());
  }

  std::size_t getRank() const
  {
    return getLeft()->getRank();
  }

  std::size_t getDimension() const
  {
    return getLeft()->getDimension();
  }
};

}

}

#endif
