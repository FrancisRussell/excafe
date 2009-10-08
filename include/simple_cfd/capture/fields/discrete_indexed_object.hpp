#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_INDEXED_OBJECT_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_INDEXED_OBJECT_HPP

#include <memory>
#include "discrete_traits.hpp"
#include "indexable_value.hpp"
#include "discrete_expr_visitor.hpp"
#include "temporal_index_expr.hpp"
#include <simple_cfd/exception.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <simple_cfd/capture/indices/propagation_rule.hpp>
#include <simple_cfd/capture/indices/index_propagation_all.hpp>
#include <simple_cfd/capture/indices/index_propagation_except.hpp>

namespace cfd
{

namespace detail
{

template<typename discrete_object_tag>
class AbstractDiscreteObjectIndexed : public DiscreteTraits<discrete_object_tag>::expr_t
{
protected:
  typedef IndexableValue<discrete_object_tag> parent_t; 
  typedef typename parent_t::value_ptr parent_ptr;
  parent_ptr parent;
  TemporalIndexExpr indexExpr;

  class IndexValidator : public boost::static_visitor<void>
  {
  private:
    const signed offsetValue;

  public:
    IndexValidator(const signed _offsetValue) : offsetValue(_offsetValue)
    {
    }

    void operator()(const TemporalIndexOffset::absolute_tag&) const
    {
      CFD_EXCEPTION("Cannot use an absolute value when indexing an r-value");
    }

    void operator()(const TemporalIndexOffset::relative_tag&) const
    {
      if(offsetValue > 0)
        CFD_EXCEPTION("Cannot refer to values from a future iteration");
    }

    void operator()(const TemporalIndexOffset::final_tag&) const
    {
      if (offsetValue > 0)
        CFD_EXCEPTION("Cannot refer to values beyond the final iteration of a loop");
    }
  };

  class InsideLoopHelper : public boost::static_visitor<bool>
  {
  public:
    bool operator()(const TemporalIndexOffset::absolute_tag&) const
    {
      return false;
    }

    bool operator()(const TemporalIndexOffset::relative_tag&) const
    {
      return true;
    }

    bool operator()(const TemporalIndexOffset::final_tag&) const
    {
      return false;
    }
  };

  std::auto_ptr<PropagationRule> constructRule(DiscreteExpr& expr, const bool insideLoop)
  {
    if (insideLoop)
      return std::auto_ptr<PropagationRule>(new IndexPropagationAll(expr, *this));
    else
      return std::auto_ptr<PropagationRule>(new IndexPropagationExcept(expr, *this, *parent->getIndexVariable()));
  }

  bool isInsideLoop() const
  {
    TemporalIndexOffset offset = indexExpr.getOffset();
    TemporalIndexOffset::offset_t offsetType = offset.getType();

    const InsideLoopHelper helper;
    return boost::apply_visitor(helper, offsetType);
  }

  signed getOffsetValue() const
  {
    const TemporalIndexOffset offset = indexExpr.getOffset();
    return offset.getValue();
  }

public:
  AbstractDiscreteObjectIndexed(const parent_ptr& _parent, const TemporalIndexExpr& _indexExpr) :
    parent(_parent), indexExpr(_indexExpr)
  {
    if(indexExpr.getIndex() != parent->getIndexVariable())
      CFD_EXCEPTION("Incorrect index used on expression rhs");
    
    TemporalIndexOffset offset = indexExpr.getOffset();
    TemporalIndexOffset::offset_t offsetType = offset.getType();

    const IndexValidator validator(offset.getValue());
    boost::apply_visitor(validator, offsetType);
  }

  TemporalIndexExpr getTemporalIndexExpr() const
  {
    return indexExpr;
  }

  parent_t& getParent() const
  {
    return *parent;
  }

  TemporalIndexSet getTemporalIndices() const
  {
    TemporalIndexSet indices;

    if (isInsideLoop())
      indices += &(*parent->getIndexVariable());

    return indices;
  }

  TemporalIndexSet getLoopDependencies() const
  {
    TemporalIndexSet indices;

    // Note: !isInsideLoop -> final
    // This class should never exist for an absolute reference
    if (!isInsideLoop())
      indices += &(*parent->getIndexVariable());

    return indices;
  }

  virtual PropagationRules getPropagationRules()
  {
    const bool insideLoop = isInsideLoop();
    PropagationRules rules;
    rules.insert(constructRule(*parent->getIterationAssignment(), insideLoop));

    for(typename IndexableValue<discrete_object_tag>::init_iterator initIter(parent->begin_inits());
      initIter!=parent->end_inits(); ++initIter)
    {
      rules.insert(constructRule(*initIter, insideLoop));
    }

    return rules;
  }

  virtual std::set<DiscreteExpr*> getDependencies() const 
  {
    std::set<DiscreteExpr*> dependencies;

    if (isInsideLoop() && getOffsetValue() == 0)
      dependencies.insert(&(*parent->getIterationAssignment()));

    return dependencies;
  }
};

// FIXME: replicating the register and unregister calls in each subclass is dangerous.

class DiscreteIndexedScalar : public AbstractDiscreteObjectIndexed<discrete_scalar_tag>
{
public:
  DiscreteIndexedScalar(const parent_ptr& _parent, const TemporalIndexExpr& indexExpr) :
    AbstractDiscreteObjectIndexed<discrete_scalar_tag>(_parent, indexExpr)
  {
    this->parent->registerReference(*this);
  }

  virtual void accept(DiscreteExprVisitor& visitor)
  {
    visitor.visit(*this);
  }

  ~DiscreteIndexedScalar()
  {
    this->parent->unregisterReference(*this);
  }
};

class DiscreteIndexedField : public AbstractDiscreteObjectIndexed<discrete_field_tag>
{
private:
  typedef DiscreteTraits<discrete_field_tag>::expr_ptr expr_ptr;

public:
  DiscreteIndexedField(const parent_ptr& _parent, const TemporalIndexExpr& indexExpr) :
    AbstractDiscreteObjectIndexed<discrete_field_tag>(_parent, indexExpr)
  {
    this->parent->registerReference(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    const expr_ptr iterationAssignment = this->parent->getIterationAssignment();
    return iterationAssignment->getFunctionSpace();
  }

  virtual void accept(DiscreteExprVisitor& visitor)
  {
    visitor.visit(*this);
  }

  ~DiscreteIndexedField()
  {
    this->parent->unregisterReference(*this);
  }
};

class DiscreteIndexedOperator : public AbstractDiscreteObjectIndexed<discrete_operator_tag>
{
public:
  DiscreteIndexedOperator(const parent_ptr& _parent, const TemporalIndexExpr& indexExpr) :
    AbstractDiscreteObjectIndexed<discrete_operator_tag>(_parent, indexExpr)
  {
    this->parent->registerReference(*this);
  }

  virtual void accept(DiscreteExprVisitor& visitor)
  {
    visitor.visit(*this);
  }

  ~DiscreteIndexedOperator()
  {
    this->parent->unregisterReference(*this);
  }
};

}

}

#endif
