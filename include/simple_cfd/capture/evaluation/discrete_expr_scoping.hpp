#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_HPP

#include <map>
#include <set>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <boost/variant.hpp>
#include <boost/variant/apply_visitor.hpp>
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
  typedef boost::variant<DiscreteExpr*, TemporalIndexValue*> evaluatable_t;

  std::map<TemporalIndexValue*, DiscreteExprScoping> loops;
  std::set<DiscreteExpr*> exprs;

  // TODO: remove replicated code from operators
  class OrderCalculationHelper : public boost::static_visitor<void>
  {
  private:
    DiscreteExprScoping& parent;
    std::set<evaluatable_t>& visited;
    std::vector<evaluatable_t>& ordered;

  public:
    OrderCalculationHelper(DiscreteExprScoping& _parent, std::set<evaluatable_t>& _visited,
      std::vector<evaluatable_t>& _ordered) : parent(_parent), visited(_visited), ordered(_ordered)
    {
    }

    void operator()(DiscreteExpr* const expr) const
    {
      if (visited.find(expr) != visited.end()) return;

      const std::set<DiscreteExpr*> dependencies = expr->getDependencies();
      for(std::set<DiscreteExpr*>::const_iterator depIter(dependencies.begin()); depIter!=dependencies.end(); ++depIter)
      {
        evaluatable_t evaluatable(*depIter);
        boost::apply_visitor(*this, evaluatable);
      }

      const TemporalIndexSet loopDependencies = expr->getLoopDependencies();
      assert(loopDependencies.size() == 0 || loopDependencies.size() == 1);

      for(TemporalIndexSet::const_iterator loopIter(loopDependencies.begin()); loopIter!=loopDependencies.end(); ++loopIter)
      {
        evaluatable_t evaluatable(&(*loopIter));
        boost::apply_visitor(*this, evaluatable);
      }

      ordered.push_back(expr);
      const bool inserted = visited.insert(expr).second;

      if (!inserted)
        CFD_EXCEPTION("Cycle detected in discrete expression DAG.");
    }

    void operator()(TemporalIndexValue* loopIndex) const
    {
      if (visited.find(loopIndex) != visited.end()) return;

      const std::map<TemporalIndexValue*, DiscreteExprScoping>::iterator loopIter = parent.loops.find(loopIndex);
      assert(loopIter != parent.loops.end());

      DiscreteExprScoping& loopScope = loopIter->second;
      const std::set<DiscreteExpr*> dependencies = loopScope.getDependencies();

      for(std::set<DiscreteExpr*>::const_iterator depIter(dependencies.begin()); depIter!=dependencies.end(); ++depIter)
      {
        evaluatable_t evaluatable(*depIter);
        boost::apply_visitor(*this, evaluatable);
      }

      ordered.push_back(loopIndex);
      const bool inserted = visited.insert(loopIndex).second;

      if (!inserted)
        CFD_EXCEPTION("Cycle detected in discrete expression DAG.");
    }
  };

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

  std::set<DiscreteExpr*> getDependencies() const
  {
    std::set<DiscreteExpr*> dependencies;

    for(std::set<DiscreteExpr*>::const_iterator exprIter(exprs.begin()); exprIter!=exprs.end(); ++exprIter)
    {
      const std::set<DiscreteExpr*> exprDependencies((*exprIter)->getDependencies());
      dependencies.insert(exprDependencies.begin(), exprDependencies.end());
    }

    for (std::map<TemporalIndexValue*, DiscreteExprScoping>::const_iterator loopIter(loops.begin());
      loopIter!=loops.end(); ++loopIter)
    {
      const std::set<DiscreteExpr*> loopDependencies(loopIter->second.getDependencies());
      dependencies.insert(loopDependencies.begin(), loopDependencies.end());
    }

    // Now we've found all dependencies of the nodes and loops in this scope, we need to subtract
    // any nodes actually present in this scope.

    std::set<DiscreteExpr*> subtractedDependencies;
    std::set_difference(dependencies.begin(), dependencies.end(), exprs.begin(), exprs.end(), 
      std::inserter(subtractedDependencies, subtractedDependencies.begin()));

    return subtractedDependencies;
  }

  void order(const std::set<DiscreteExpr*>& wanted)
  {
    std::set<evaluatable_t> visited;
    std::vector<evaluatable_t> ordered;
    OrderCalculationHelper helper(*this, visited, ordered);

    for(std::set<DiscreteExpr*>::const_iterator wantedIter(wanted.begin()); wantedIter!=wanted.end(); ++wantedIter)
    {
      evaluatable_t evaluatable(*wantedIter);
      boost::apply_visitor(helper, evaluatable);
    }

    //TODO: We also need to order all nested loops
  }
};

}

}
#endif
