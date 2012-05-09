#include <map>
#include <set>
#include <cstddef>
#include <algorithm>
#include <utility>
#include <iterator>
#include <boost/variant/apply_visitor.hpp>
#include <excafe/exception.hpp>
#include <excafe/capture/capture_fwd.hpp>
#include <excafe/capture/fields/temporal_index_value.hpp>
#include <excafe/capture/fields/indexable_value.hpp>
#include <excafe/capture/fields/temporal_index_set.hpp>
#include <excafe/capture/fields/discrete_expr.hpp>
#include <excafe/capture/evaluation/evaluation_fwd.hpp>
#include <excafe/capture/evaluation/expression_values.hpp>
#include <excafe/capture/evaluation/discrete_expr_scoping.hpp>
#include <excafe/capture/evaluation/discrete_expr_scoping_visitor.hpp>

namespace excafe
{

namespace detail
{

TemporalIndexSet DiscreteExprScoping::getLoopIndices() const
{
  TemporalIndexSet indices;

  for(std::map<TemporalIndexValue*, DiscreteExprScoping>::const_iterator loopIter(loops.begin());
    loopIter!=loops.end(); ++loopIter)
  {
    indices += loopIter->first;
  }

  return indices;
}

std::size_t DiscreteExprScoping::maxIndices(const std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices)
{
  std::size_t max = 0;
  for(std::map<DiscreteExpr*, TemporalIndexSet>::const_iterator exprIter(exprIndices.begin());
    exprIter!=exprIndices.end(); ++exprIter)
  {
    max = std::max(max, exprIter->second.size());
  }
  return max;
}

DiscreteExprScoping& DiscreteExprScoping::getLoop(TemporalIndexValue& index)
{
  const std::map<TemporalIndexValue*, DiscreteExprScoping>::iterator loopIter(loops.find(&index));

  if (loopIter!=loops.end())
  {
    return loopIter->second;
  }
  else
  {
    return loops.insert(std::make_pair(&index, DiscreteExprScoping(index))).first->second;
  }
}

/* This function allows expressions to be added to this scope. To be able to infer the loop nesting
   structure, expressions must be added in order of increasing numbers of indices.
*/
void DiscreteExprScoping::addExpressionNode(const TemporalIndexSet& indices, DiscreteExpr& expr)
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
    TemporalIndexSet empty;
    getLoop(index).addExpressionNode(empty, expr);
  }
  else
  {
    const TemporalIndexSet loopIndices = getLoopIndices();
    TemporalIndexSet intersection = indices.intersection(loopIndices);

    if (intersection.size() == 1)
    {
      TemporalIndexValue& index = *intersection.begin();
      const TemporalIndexSet remainingIndices = indices - &index;
      getLoop(index).addExpressionNode(remainingIndices, expr);
    }
    else
    {
      // We cannot infer the scope of this expression. No expressions were added before it
      // that enabled us to determine what loop it belongs to.
      CFD_EXCEPTION("Unable to infer loop nesting structure of discrete expression node.");
    }
  }
}

void DiscreteExprScoping::addExpressionNodes(const std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices)
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

std::set<DiscreteExpr*> DiscreteExprScoping::getDependencies() const
{
  std::set<DiscreteExpr*> dependencies;

  // Get all dependencies of expression nodes in this scope
  for(std::set<DiscreteExpr*>::const_iterator exprIter(exprs.begin()); exprIter!=exprs.end(); ++exprIter)
  {
    const std::set<DiscreteExpr*> exprDependencies((*exprIter)->getDependencies());
    dependencies.insert(exprDependencies.begin(), exprDependencies.end());
  }

  // Get all dependencies of loops in this scope
  for (std::map<TemporalIndexValue*, DiscreteExprScoping>::const_iterator loopIter(loops.begin());
    loopIter!=loops.end(); ++loopIter)
  {
    const std::set<DiscreteExpr*> loopDependencies(loopIter->second.getDependencies());
    dependencies.insert(loopDependencies.begin(), loopDependencies.end());
  }

  // Get all dependencies of initialisers required by this loop
  if (!isGlobalScope())
  {
    addInitialisers(thisLoopIndex->getIndexableScalars(), dependencies);
    addInitialisers(thisLoopIndex->getIndexableFields(), dependencies);
    addInitialisers(thisLoopIndex->getIndexableOperators(), dependencies);
  }

  // Now we've found all dependencies of the nodes and loops in this scope, we need to subtract
  // any nodes actually present in this scope.

  std::set<DiscreteExpr*> subtractedDependencies;
  std::set_difference(dependencies.begin(), dependencies.end(), exprs.begin(), exprs.end(), 
    std::inserter(subtractedDependencies, subtractedDependencies.begin()));

  return subtractedDependencies;
}

void DiscreteExprScoping::order(const std::set<DiscreteExpr*>& wanted)
{
  OrderCalculationHelper helper(*this);

  // Determine dependencies for wanted values
  for(std::set<DiscreteExpr*>::const_iterator wantedIter(wanted.begin()); wantedIter!=wanted.end(); ++wantedIter)
  {
    evaluatable_t evaluatable(*wantedIter);
    boost::apply_visitor(helper, evaluatable);
  }

  if (!isGlobalScope())
  {
    // Determine dependencies for iteration assignment of indexable values
    orderIndexable(thisLoopIndex->getIndexableScalars(), helper);
    orderIndexable(thisLoopIndex->getIndexableFields(), helper);
    orderIndexable(thisLoopIndex->getIndexableOperators(), helper);

    // Determine dependencies for termination condition
    evaluatable_t evaluatable(&thisLoopIndex->getTermination());
    boost::apply_visitor(helper, evaluatable);
  }

  ordered = helper.getOrdered();

  // Now order nested loops
  const std::set<DiscreteExpr*> empty;
  for(std::map<TemporalIndexValue*, DiscreteExprScoping>::iterator loopIter(loops.begin());
    loopIter!=loops.end(); ++loopIter)
  {
    loopIter->second.order(empty);
  }
}

void DiscreteExprScoping::accept(DiscreteExprScopingVisitor& v)
{
  if (isGlobalScope())
  {
    v.visitBlock(*this);
  }
  else
  {
    v.visitLoop(*this, thisLoopIndex);
  }
}

boost::shared_ptr<TemporalIndexValue> DiscreteExprScoping::globalScope(new TemporalIndexValue());

}

}
