#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_EVALUATION_STRATEGY_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_EVALUATION_STRATEGY_HPP

#include <map>
#include <set>
#include <simple_cfd/capture/fields/fields_fwd.hpp>
#include <simple_cfd/capture/fields/discrete_expr_container.hpp>
#include <simple_cfd/capture/fields/discrete_expr_set.hpp>
#include <simple_cfd/capture/fields/temporal_index_set.hpp>
#include "expression_values.hpp"
#include "evaluation_visitor.hpp"
#include "discrete_expr_scoping.hpp"

namespace cfd
{

namespace detail
{

class EvaluationStrategy
{
private:
  const std::set<DiscreteExpr*> wantedExprs;
  const DiscreteExprContainer expr;
  DiscreteExprScoping scoping;

  template<typename discrete_object_tag>
  static void sortExpressions(std::map<TemporalIndexSet, DiscreteExprContainer>& indicesToExprMap,
    const DiscreteExprSet<discrete_object_tag>& expressions)
  {
    for(typename DiscreteExprSet<discrete_object_tag>::expr_iter exprIter(expressions.begin_expr()); 
      exprIter!=expressions.end_expr(); ++exprIter)
    {
      indicesToExprMap[exprIter->getTemporalIndices()].insert(*exprIter);
    }
  }

  template<typename discrete_object_tag>
  static void addIndices(std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices, const DiscreteExprSet<discrete_object_tag>& expressions)
  {
    for(typename DiscreteExprSet<discrete_object_tag>::expr_iter exprIter(expressions.begin_expr()); 
      exprIter!=expressions.end_expr(); ++exprIter)
    {
      exprIndices[&(*exprIter)] += exprIter->getTemporalIndices();
    }
  }

  std::map<DiscreteExpr*, TemporalIndexSet> findExpressionIndices() const;

public:
  EvaluationStrategy(const DiscreteExprContainer& _expr, const std::set<DiscreteExpr*>& _wantedExprs);
  void buildExprScoping();

  template<std::size_t D>
  ExpressionValues<D> execute(Scenario<D>& scenario)
  {
    EvaluationVisitor<D> evaluationVisitor(scenario);
    scoping.accept(evaluationVisitor);
    return evaluationVisitor.getWanted();
  }
};

}

}

#endif
