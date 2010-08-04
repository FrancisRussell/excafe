#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_HPP

#include <set>
#include "fields_fwd.hpp"
#include <simple_cfd/capture/indices/indices_fwd.hpp>

namespace cfd
{

namespace detail
{

class DiscreteExpr
{
public:
  virtual void accept(DiscreteExprVisitor& visitor) = 0;
  virtual TemporalIndexSet getTemporalIndices() const;
  virtual TemporalIndexSet getLoopDependencies() const;
  virtual std::set<DiscreteExpr*> getDependencies() const = 0;
  virtual PropagationRules getPropagationRules() = 0;
  virtual ~DiscreteExpr() {}
};

}

}

#endif