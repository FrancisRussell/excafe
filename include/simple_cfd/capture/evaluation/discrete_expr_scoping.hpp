#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_HPP

#include <map>
#include <set>
#include <cstddef>
#include <algorithm>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/capture/fields/temporal_index_set.hpp>
#include <simple_cfd/capture/fields/discrete_expr.hpp>


namespace cfd
{

namespace detail
{

class DiscreteExprScoping
{
private:
  std::map<TemporalIndexValue*, DiscreteExprScoping> loops;
  std::set<DiscreteExpr*> exprs;

  TemporalIndexSet getLoopIndices() const
  {
    TemporalIndexSet indices;

    for(std::map<TemporalIndexValue*, DiscreteExprScoping>::const_iterator loopIter(loops.begin());
      loopIter!=loops.end(); ++loopIter)
    {
      indices += loopIter->first;
    }

    return indices;
  }

  static std::size_t maxIndices(const std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices)
  {
    std::size_t max = 0;
    for(std::map<DiscreteExpr*, TemporalIndexSet>::const_iterator exprIter(exprIndices.begin());
      exprIter!=exprIndices.end(); ++exprIter)
    {
      max = std::max(max, exprIter->second.size());
    }
    return max;
  }

public:
  /* This function allows expressions to be added to this scope. To be able to infer the loop nesting
     structure, expressions must be added in order of increasing numbers of indices.
  */
  void addExpressionNode(const TemporalIndexSet& indices, DiscreteExpr& expr)
  {
    if (indices.size() == 0)
    {
      // The expression belongs in this scope
      exprs.insert(&expr);
    }
    else if (indices.size() == 1)
    {
      // The expression belongs in an inferable immediate sub-scope (which we MAY have to create)
      TemporalIndexValue& index = *indices.begin();
      const TemporalIndexSet empty;
      loops[&index].addExpressionNode(empty, expr);
    }
    else
    {
      const TemporalIndexSet loopIndices = getLoopIndices();
      TemporalIndexSet intersection = indices.intersection(loopIndices);

      if (intersection.size() == 1)
      {
        TemporalIndexValue& index = *intersection.begin();
        const TemporalIndexSet remainingIndices = indices - &index;
        loops[&index].addExpressionNode(remainingIndices, expr);
      }
      else
      {
        // We cannot infer the scope of this expression. No expressions were added before it
        // that enabled us to determine what loop it belongs to.
        CFD_EXCEPTION("Unable to infer loop nesting structure of discrete expression node.");
      }
    }
  }

  void addExpressionNodes(const std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices)
  {
    const std::size_t numIndices = maxIndices(exprIndices);

    for(std::size_t indexDepth = 0; indexDepth <=numIndices; ++indexDepth)
    {
      for(std::map<DiscreteExpr*, TemporalIndexSet>::const_iterator exprIter(exprIndices.begin());
        exprIter!=exprIndices.end(); ++exprIter)
      {
        if (exprIter->second.size() == indexDepth)
        {
          addExpressionNode(exprIter->second, *exprIter->first);
        }
      }
    }
  }
};

}

}
#endif
