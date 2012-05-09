#ifndef EXCAFE_CAPTURE_FIELDS_DISCRETE_EXPR_HPP
#define EXCAFE_CAPTURE_FIELDS_DISCRETE_EXPR_HPP

#include <set>
#include "fields_fwd.hpp"
#include <excafe/capture/indices/indices_fwd.hpp>

namespace excafe
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
