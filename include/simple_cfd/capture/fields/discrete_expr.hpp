#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_HPP

namespace cfd
{

namespace detail
{

class DiscreteExpr
{
public:
  virtual void accept(DiscreteExprVisitor& visitor) = 0;
  virtual TemporalIndexSet getTemporalIndices() const = 0;
  virtual ~DiscreteExpr() {}
};

}

}

#endif
