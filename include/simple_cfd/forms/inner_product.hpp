#ifndef SIMPLE_CFD_FORMS_INNER_PRODUCT_HPP
#define SIMPLE_CDS_FORMS_INNER_PRODUCT_HPP

#include <cstddef>
#include "binary_operator.hpp"

namespace cfd
{

namespace forms
{

class InnerProduct : public BinaryOperator
{
public:
  InnerProduct(Field::reference_t l, Field::reference_t r) : BinaryOperator(l, r)
  {
    assert(l->getRank() + r->getRank() >= 2);
    assert(l->getDimension() == r->getDimension());
  }

  std::size_t getRank() const
  {
    return getLeft()->getRank() + getRight()->getRank() - 2;
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
