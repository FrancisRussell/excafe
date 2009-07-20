#ifndef SIMPLE_CFD_FORMS_OUTER_PRODUCT_HPP
#define SIMPLE_CFD_FORMS_OUTER_PRODUCT_HPP

#include <cstddef>
#include "binary_operator.hpp"

namespace cfd
{

namespace forms
{

class OuterProduct : public BinaryOperator
{
public:
  OuterProduct(Field::reference_t l, Field::reference_t r) : BinaryOperator(l, r)
  {
    assert(l->getDimension() == r->getDimension());
  }

  std::size_t getRank() const
  {
    return getLeft()->getRank() + getRight()->getRank();
  }

  std::size_t getDimension() const
  {
    return getLeft()->getDimension();
  }
  
  virtual void accept(FieldVisitor& v)
  {
    v.enter(*this);
    getLeft()->accept(v);
    getRight()->accept(v);
    v.exit(*this);
  }
};

}

}

#endif
