#ifndef SIMPLE_CFD_FORMS_FACET_NORMAL_HPP
#define SIMPLE_CFD_FORMS_FACET_NORMAL_HPP

#include <cstddef>
#include "field.hpp"
#include "field_visitor.hpp"

namespace cfd
{

namespace forms
{

class facet_normal_tag {};

class FacetNormal : public Field
{
public:
  FacetNormal()
  {
  }

  std::size_t getRank() const
  {
    return 1;
  }

  virtual void accept(FieldVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
