#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_HPP

#include <map>
#include <set>
#include <vector>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <boost/variant.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/capture/fields/temporal_index_value.hpp>
#include <simple_cfd/capture/fields/indexable_value.hpp>
#include <simple_cfd/capture/fields/temporal_index_set.hpp>
#include <simple_cfd/capture/fields/discrete_expr.hpp>

namespace cfd
{

namespace detail
{

class DiscreteExprScoping
{
private:
  static boost::shared_ptr<TemporalIndexValue> globalScope;
  typedef boost::variant<DiscreteExpr*, TemporalIndexValue*> evaluatable_t;

  TemporalIndexValue* const thisLoopIndex;
  std::map<TemporalIndexValue*, DiscreteExprScoping> loops;
  std::set<DiscreteExpr*> exprs;

  DiscreteExprScoping(TemporalIndexValue& _thisLoopIndex) : thisLoopIndex(&_thisLoopIndex)
  {
  }

  DiscreteExprScoping& getLoop(TemporalIndexValue& index)
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

  bool isGlobalScope() const
  {
    return thisLoopIndex == globalScope.get();
  }

  template<typename discrete_object_tag>
  void orderIndexable(const std::set<IndexableValue<discrete_object_tag>*>& indexableValues, OrderCalculationHelper& helper)
  {
    typedef typename std::set<IndexableValue<discrete_object_tag>*>::const_iterator indexable_set_iter;
    for(indexable_set_iter indexableIter(indexableValues.begin()); indexableIter!=indexableValues.end(); ++indexableIter)
    {
      evaluatable_t evaluatable(&(*(*indexableIter)->getIterationAssignment()));
      boost::apply_visitor(helper, evaluatable);
    }
  }

  
  template<typename discrete_object_tag>
  void addInitialisers(const std::set<IndexableValue<discrete_object_tag>*>& indexableValues, std::set<DiscreteExpr*> initialisers) const
  {
    typedef typename std::set<IndexableValue<discrete_object_tag>*>::const_iterator indexable_set_iter;
    typedef typename IndexableValue<discrete_object_tag>::init_iterator indexable_init_iter;

    for(indexable_set_iter indexableIter(indexableValues.begin()); indexableIter!=indexableValues.end(); ++indexableIter)
    {
      for(indexable_init_iter initIter((*indexableIter)->begin_inits()); initIter!=(*indexableIter)->end_inits(); ++initIter)
      {
        initialisers.insert(&(*initIter));
      }
    }
  }

public:
  DiscreteExprScoping() : thisLoopIndex(globalScope.get())
  {
  }

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

  void order(const std::set<DiscreteExpr*>& wanted)
  {
    std::set<evaluatable_t> visited;
    std::vector<evaluatable_t> ordered;
    OrderCalculationHelper helper(*this, visited, ordered);

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

    // Now order nested loops
    const std::set<DiscreteExpr*> empty;
    for(std::map<TemporalIndexValue*, DiscreteExprScoping>::iterator loopIter(loops.begin());
      loopIter!=loops.end(); ++loopIter)
    {
      loopIter->second.order(empty);
    }
  }
};

}

}
#endif
