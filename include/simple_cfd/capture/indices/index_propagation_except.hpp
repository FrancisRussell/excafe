#ifndef SIMPLE_CFD_CAPTURE_INDICES_INDEX_PROPAGATION_EXCEPT_HPP
#define SIMPLE_CFD_CAPTURE_INDICES_INDEX_PROPAGATION_EXCEPT_HPP

#include "propagation_rule_visitor.hpp"

namespace cfd
{

namespace detail
{

class IndexPropagationExcept : public PropagationRule
{
private:
  TemporalIndexValue* const index;

public:
  IndexPropagationExcept(DiscreteExpr& from, DiscreteExpr& to, TemporalIndexValue& i) : 
    PropagationRule(from, to), index(&i)
  {
  }

  TemporalIndexValue& getExcludedIndex() const
  {
    return *index;
  }

  virtual void accept(PropagationRuleVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
