#ifndef EXCAFE_FORMS_FACET_NORMAL_HPP
#define EXCAFE_FORMS_FACET_NORMAL_HPP

#include <cstddef>
#include "field_expr.hpp"
#include "field_visitor.hpp"

namespace excafe
{

namespace detail
{

class FacetNormal : public FieldExpr
{
public:
  FacetNormal()
  {
  }

  virtual void accept(FieldVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
