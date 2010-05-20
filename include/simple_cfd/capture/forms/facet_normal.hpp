#ifndef SIMPLE_CFD_FORMS_FACET_NORMAL_HPP
#define SIMPLE_CFD_FORMS_FACET_NORMAL_HPP

#include <cstddef>
#include "field_expr.hpp"
#include "field_visitor.hpp"

namespace cfd
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
