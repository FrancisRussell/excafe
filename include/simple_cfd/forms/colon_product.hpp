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
    assert(l->getRank() == r->getRank());
  }

  std::size_t getRank() const
  {
    return 0;
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
