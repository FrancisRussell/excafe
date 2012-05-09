#ifndef EXCAFE_CAPTURE_FIELDS_SCALAR_EXPR_HPP
#define EXCAFE_CAPTURE_FIELDS_SCALAR_EXPR_HPP

#include "fields_fwd.hpp"
#include "discrete_expr.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

namespace excafe
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
