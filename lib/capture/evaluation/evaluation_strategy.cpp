#include <map>
#include <set>
#include <simple_cfd/capture/fields/fields_fwd.hpp>
#include <simple_cfd/capture/fields/temporal_index_set.hpp>
#include <simple_cfd/capture/evaluation/evaluation_strategy.hpp>
#include <simple_cfd/capture/fields/discrete_expr_container.hpp>

namespace cfd
{

namespace detail
{

std::map<DiscreteExpr*, TemporalIndexSet> EvaluationStrategy::findExpressionIndices() const
{
  std::map<DiscreteExpr*, TemporalIndexSet> exprIndices;
  addIndices(exprIndices, expr.getScalarExpressions());
  addIndices(exprIndices, expr.getFieldExpressions());
  addIndices(exprIndices, expr.getOperatorExpressions());

  PropagationRules rules;
  for(std::map<DiscreteExpr*, TemporalIndexSet>::const_iterator exprIter(exprIndices.begin());
    exprIter!=exprIndices.end(); ++exprIter)
  {
    rules += exprIter->first->getPropagationRules();
  }

  rules.propagateIndices(exprIndices);
  return exprIndices;
}


EvaluationStrategy::EvaluationStrategy(const DiscreteExprContainer& _expr, const std::set<DiscreteExpr*>& _wantedExprs) :
  wantedExprs(_wantedExprs), expr(_expr)
{
  buildExprScoping();
}

void EvaluationStrategy::buildExprScoping()
{
  const std::map<DiscreteExpr*, TemporalIndexSet> exprIndices = findExpressionIndices();
  scoping.addExpressionNodes(exprIndices);
  scoping.order(wantedExprs);
}

}

}
