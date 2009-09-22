#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_EVALUATION_STRATEGY_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_EVALUATION_STRATEGY_HPP

#include <map>
#include <iostream>
#include <simple_cfd/capture/fields/discrete_expr_container.hpp>
#include <simple_cfd/capture/fields/temporal_index_set.hpp>

namespace cfd
{

namespace detail
{

class EvaluationStrategy
{
private:
  const DiscreteExprContainer expr;

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

public:
  EvaluationStrategy(const DiscreteExprContainer& _expr) : expr(_expr)
  {
    buildExprScoping();
  }

  void buildExprScoping()
  {
    std::map<TemporalIndexSet, DiscreteExprContainer> indicesToExprMap;
    sortExpressions(indicesToExprMap, expr.getScalarExpressions());

    std::cout << "Number of distinct scopes: " << indicesToExprMap.size() << std::endl;
  }
};

}

}

#endif
