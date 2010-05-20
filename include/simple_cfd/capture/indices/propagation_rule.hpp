#ifndef SIMPLE_CFD_CAPTURE_INDICES_PROPAGATION_RULE_HPP
#define SIMPLE_CFD_CAPTURE_INDICES_PROPAGATION_RULE_HPP

#include "indices_fwd.hpp"
#include <simple_cfd/capture/fields/discrete_expr.hpp>

namespace cfd
{

namespace detail
{

class PropagationRule
{
private:
  DiscreteExpr* const from;
  DiscreteExpr* const to;

public:
  PropagationRule(DiscreteExpr& _from, DiscreteExpr& _to) : from(&_from), to(&_to)
  {
  }

  DiscreteExpr& getFrom() const
  {
    return *from;
  }

  DiscreteExpr& getTo() const
  {
    return *to;
  }

  virtual void accept(PropagationRuleVisitor& v) = 0;

  virtual ~PropagationRule() {}
};

}

}

#endif
