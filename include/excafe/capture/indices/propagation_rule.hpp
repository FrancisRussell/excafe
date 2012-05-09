#ifndef EXCAFE_CAPTURE_INDICES_PROPAGATION_RULE_HPP
#define EXCAFE_CAPTURE_INDICES_PROPAGATION_RULE_HPP

#include "indices_fwd.hpp"
#include <excafe/capture/fields/discrete_expr.hpp>

namespace excafe
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
