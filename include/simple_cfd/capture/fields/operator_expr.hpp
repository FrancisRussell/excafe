#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_EXPR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_EXPR_HPP

#include "fields_fwd.hpp"
#include "function_space_expr.hpp"
#include "discrete_expr.hpp"
#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace detail
{

class OperatorExpr : public DiscreteExpr
{
public:
  typedef boost::shared_ptr<OperatorExpr> expr_ptr;

  virtual FunctionSpaceExpr::expr_ptr getTrialSpace() const = 0;
  virtual FunctionSpaceExpr::expr_ptr getTestSpace() const = 0;
};

}

}

#endif
