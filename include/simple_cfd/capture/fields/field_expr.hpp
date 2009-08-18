#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELD_EXPR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELD_EXPR_HPP

#include "fields_fwd.hpp"

namespace cfd
{

namespace detail
{

class FieldExpr
{
public:
  typedef boost::shared_ptr<FieldExpr> expr_ptr;

  virtual void accept(FieldVisitor& f) = 0;
  virtual ~FieldExpr() {}
};

}

}

#endif
