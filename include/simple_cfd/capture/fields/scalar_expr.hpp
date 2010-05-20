#ifndef SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_EXPR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_EXPR_HPP

#include "fields_fwd.hpp"
#include "discrete_expr.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

namespace cfd
{

namespace detail
{

class ScalarExpr : public DiscreteExpr
{
public:
  typedef boost::shared_ptr<ScalarExpr> expr_ptr;
  typedef boost::weak_ptr<ScalarExpr> weak_expr_ptr;
};

}

}

#endif
