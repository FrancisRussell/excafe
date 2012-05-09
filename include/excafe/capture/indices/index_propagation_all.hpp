#ifndef EXCAFE_CAPTURE_INDICES_INDEX_PROPAGATION_ALL_HPP
#define EXCAFE_CAPTURE_INDICES_INDEX_PROPAGATION_ALL_HPP

#include "propagation_rule.hpp"
#include "propagation_rule_visitor.hpp"

namespace excafe
{

namespace detail
{

class IndexPropagationAll : public PropagationRule
{
public:
  IndexPropagationAll(DiscreteExpr& from, DiscreteExpr& to) : PropagationRule(from, to)
  {
  }

  virtual void accept(PropagationRuleVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
