#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELD_EMPTY_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELD_EMPTY_HPP

#include "field_visitor.hpp"

namespace cfd
{

namespace detail
{

class FieldEmpty : public FieldExpr
{
  void accept(FieldVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
