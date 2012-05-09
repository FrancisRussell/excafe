#ifndef EXCAFE_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_HPP
#define EXCAFE_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_HPP

#include <map>
#include <set>
#include <vector>
#include <cstddef>
#include <cassert>
#include <boost/variant.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <excafe/exception.hpp>
#include <excafe/capture/capture_fwd.hpp>
#include <excafe/capture/fields/temporal_index_value.hpp>
#include <excafe/capture/fields/indexable_value.hpp>
#include <excafe/capture/fields/temporal_index_set.hpp>
#include <excafe/capture/fields/discrete_expr.hpp>
#include <excafe/capture/evaluation/evaluation_fwd.hpp>
#include <excafe/capture/evaluation/expression_values.hpp>

namespace excafe
{

namespace detail
{

class DiscreteExprScoping
{
public:
  typedef boost::variant<DiscreteExpr*, TemporalIndexValue*> evaluatable_t;

private:
  static boost::shared_ptr<TemporalIndexValue> globalScope;

  TemporalIndexValue* const thisLoopIndex;
  std::map<TemporalIndexValue*, DiscreteExprScoping> loops;
  std::set<DiscreteExpr*> exprs;
  std::vector<evaluatable_t> ordered;

  // TODO: remove replicated code from operators
  class OrderCalculationHelper : public boost::static_visitor<void>
  {
  private:
    DiscreteExprScoping& parent;
    std::set<evaluatable_t> visited;
    std::vector<evaluatable_t> ordered;
     
    void addEvaluatable(const evaluatable_t& e)
    {
      ordered.push_back(e);
      const bool inserted = visited.insert(e).second;

      if (!inserted)
        CFD_EXCEPTION("Cycle detected in discrete expression DAG.");
    }

  public:
    OrderCalculationHelper(DiscreteExprScoping& _parent) : parent(_parent)
    {
    }

    void operator()(DiscreteExpr* const expr)
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

      addEvaluatable(expr);
    }

    void operator()(TemporalIndexValue* const loopIndex)
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

      addEvaluatable(loopIndex);
    }

    std::vector<evaluatable_t> getOrdered() const
    {
      return ordered;
    }
  };

  DiscreteExprScoping(TemporalIndexValue& _thisLoopIndex) : thisLoopIndex(&_thisLoopIndex)
  {
  }

  TemporalIndexSet getLoopIndices() const;
  static std::size_t maxIndices(const std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices);

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
  void addInitialisers(const std::set<IndexableValue<discrete_object_tag>*>& indexableValues, std::set<DiscreteExpr*>& initialisers) const
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
  typedef std::vector<evaluatable_t>::iterator iterator;
  typedef std::vector<evaluatable_t>::const_iterator const_iterator;

  DiscreteExprScoping() : thisLoopIndex(globalScope.get())
  {
  }

  iterator begin() 
  {
    return ordered.begin();
  }

  iterator end()
  {
    return ordered.end();
  }

  const_iterator begin() const
  {
    return ordered.begin();
  }

  const_iterator end() const
  {
    return ordered.end();
  }

  DiscreteExprScoping& getLoop(TemporalIndexValue& index);
  void addExpressionNode(const TemporalIndexSet& indices, DiscreteExpr& expr);
  void addExpressionNodes(const std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices);
  std::set<DiscreteExpr*> getDependencies() const;
  void order(const std::set<DiscreteExpr*>& wanted);
  void accept(DiscreteExprScopingVisitor& v);
};

}

}
#endif
