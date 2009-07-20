#ifndef SIMPLE_CFD_FORMS_COLON_PRODUCT_HPP
#define SIMPLE_CFD_FORMS_COLON_PRODUCT_HPP

#include <cstddef>
#include "binary_operator.hpp"

namespace cfd
{

namespace forms
{

class ColonProduct : public BinaryOperator
{
public:
  ColonProduct(Field::reference_t l, Field::reference_t r) : BinaryOperator(l, r)
  {
    assert(l->getDimension() == r->getDimension());
    assert(l->getRank() == r->getRank());
  }

  std::size_t getRank() const
  {
    return 0;
  }

  std::size_t getDimension() const
  {
    return getLeft()->getDimension();
  }
};

}

}

#endif
