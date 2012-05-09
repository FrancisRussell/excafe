#ifndef EXCAFE_CAPTURE_FIELDS_OPERATOR_EXPR_HPP
#define EXCAFE_CAPTURE_FIELDS_OPERATOR_EXPR_HPP

#include "fields_fwd.hpp"
#include "function_space_expr.hpp"
#include "discrete_expr.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

namespace excafe
{

namespace detail
{

class OperatorExpr : public DiscreteExpr
{
public:
  typedef boost::shared_ptr<OperatorExpr> expr_ptr;
  typedef boost::weak_ptr<OperatorExpr> weak_expr_ptr;

  virtual FunctionSpaceExpr::expr_ptr getTrialSpace() const = 0;
  virtual FunctionSpaceExpr::expr_ptr getTestSpace() const = 0;
};

}

}

#endif
