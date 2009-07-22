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
  }

  std::size_t getRank() const
  {
    return getLeft()->getRank() + getRight()->getRank();
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
