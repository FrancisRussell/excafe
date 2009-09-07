#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELD_EXPR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELD_EXPR_HPP

#include <boost/shared_ptr.hpp>
#include "fields_fwd.hpp"
#include "function_space_expr.hpp"

namespace cfd
{

namespace detail
{

class DiscreteFieldExpr
{
public:
  typedef boost::shared_ptr<DiscreteFieldExpr> expr_ptr;

  virtual void accept(DiscreteExprVisitor& f) = 0;
  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const = 0;
  virtual ~DiscreteFieldExpr() {}
};

}

}

#endif
