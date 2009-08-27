#ifndef SIMPLE_CFD_CAPTURE_FORMS_FIELD_EXPR_HPP
#define SIMPLE_CFD_CAPTURE_FORMS_FIELD_EXPR_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "field_visitor.hpp"

namespace cfd
{

namespace detail
{

class FieldExpr
{
public:
  typedef boost::shared_ptr<FieldExpr> reference_t;
  virtual void accept(FieldVisitor& visitor) = 0;
};

}

}

#endif
