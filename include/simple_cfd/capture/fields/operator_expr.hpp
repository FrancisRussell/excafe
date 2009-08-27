#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_EXPR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_EXPR_HPP

#include "fields_fwd.hpp"
#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace detail
{

class OperatorExpr
{
public:
  typedef boost::shared_ptr<OperatorExpr> expr_ptr;

  virtual void accept(DiscreteExprVisitor& v) = 0;
  virtual ~OperatorExpr() {}
};

}

}

#endif
