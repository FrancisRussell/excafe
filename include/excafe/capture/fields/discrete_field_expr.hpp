#ifndef EXCAFE_CAPTURE_FIELDS_FIELD_EXPR_HPP
#define EXCAFE_CAPTURE_FIELDS_FIELD_EXPR_HPP

#include "fields_fwd.hpp"
#include "function_space_expr.hpp"
#include "discrete_expr.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

namespace excafe
{

namespace detail
{

class DiscreteFieldExpr : public DiscreteExpr
{
public:
  typedef boost::shared_ptr<DiscreteFieldExpr> expr_ptr;
  typedef boost::weak_ptr<DiscreteFieldExpr> weak_expr_ptr;
  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const = 0;
};

}

}

#endif
