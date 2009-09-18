#ifndef SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_EXPR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_EXPR_HPP

#include "fields_fwd.hpp"
#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace detail
{

class ScalarExpr
{
public:
  typedef boost::shared_ptr<ScalarExpr> expr_ptr;
  virtual TemporalIndexSet getTemporalIndices() const = 0;
  virtual void accept(DiscreteExprVisitor& visitor) = 0;
  virtual ~ScalarExpr() {}
};

}

}

#endif
