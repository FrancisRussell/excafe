#ifndef EXCAFE_CAPTURE_FORMS_FIELD_BASIS_HPP
#define EXCAFE_CAPTURE_FORMS_FIELD_BASIS_HPP

#include <cstddef>
#include "field_expr.hpp"
#include <excafe/capture/fields/element.hpp>

namespace excafe
{

namespace detail
{

class FieldBasis : public FieldExpr
{
private:
  Element element;

public:
  FieldBasis(const Element _element) : element(_element)
  {
  }

  virtual void accept(FieldVisitor& v)
  {
    v.visit(*this);
  }

  Element getElement() const
  {
    return element;
  }
};

}

}
#endif
