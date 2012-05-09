#ifndef EXCAFE_CAPTURE_FORMS_FIELD_EXPR_HPP
#define EXCAFE_CAPTURE_FORMS_FIELD_EXPR_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "field_visitor.hpp"

namespace excafe
{

namespace detail
{

class FieldExpr
{
public:
  typedef boost::shared_ptr<FieldExpr> reference_t;
  virtual void accept(FieldVisitor& visitor) = 0;
  virtual ~FieldExpr() {}
};

}

}

#endif
