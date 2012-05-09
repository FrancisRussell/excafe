#include <map>
#include <excafe/capture/fields/temporal_index_set.hpp>
#include <excafe/capture/indices/index_propagation_all.hpp>
#include <excafe/capture/indices/index_propagation_except.hpp>
#include <excafe/capture/indices/propagation_rules.hpp>

namespace excafe
{

namespace detail
{

void RulePropagationHelper::visit(IndexPropagationAll& i)
{
  const std::map<DiscreteExpr*, TemporalIndexSet>::const_iterator fromIter(exprIndices.find(&i.getFrom()));
  const std::map<DiscreteExpr*, TemporalIndexSet>::iterator toIter(exprIndices.find(&i.getTo()));

  assert(fromIter != exprIndices.end());
  assert(toIter != exprIndices.end());

  toIter->second += fromIter->second;
}

void RulePropagationHelper::visit(IndexPropagationExcept& i)
{
  const std::map<DiscreteExpr*, TemporalIndexSet>::const_iterator fromIter(exprIndices.find(&i.getFrom()));
  const std::map<DiscreteExpr*, TemporalIndexSet>::iterator toIter(exprIndices.find(&i.getTo()));

  assert(fromIter != exprIndices.end());
  assert(toIter != exprIndices.end());

  TemporalIndexSet indices(fromIter->second);
  indices -= &i.getExcludedIndex();
  toIter->second += indices;
}

void PropagationRules::propagateIndices(std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices) const
{
  bool converged = false;
  RulePropagationHelper propagationHelper(exprIndices);

  while(!converged)
  {
    const std::map<DiscreteExpr*, TemporalIndexSet> exprIndicesOld(exprIndices);

    for(const_iterator ruleIter(begin()); ruleIter!=end(); ++ruleIter)
    {
      ruleIter->accept(propagationHelper);
    }

    converged = (exprIndices == exprIndicesOld);
  }
}

PropagationRules& PropagationRules::operator+=(const PropagationRules& p)
{
  rules.insert(p.rules.begin(), p.rules.end());
  return *this;
}

}
 
}
