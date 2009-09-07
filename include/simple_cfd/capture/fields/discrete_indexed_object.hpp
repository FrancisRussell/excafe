#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_OBJECT_INDEXED_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_OBJECT_INDEXED_HPP

#include <boost/static_assert.hpp>
#include "discrete_traits.hpp"
#include "indexable_value.hpp"
#include "discrete_expr_visitor.hpp"

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
  TemporalIndexExpr index;

public:
  AbstractDiscreteObjectIndexed(const parent_ptr& _parent, const TemporalIndexExpr indexExpr) :
    parent(_parent), index(indexExpr)
  {
  }
};

class DiscreteIndexedScalar : public AbstractDiscreteObjectIndexed<discrete_scalar_tag>
{
public:
  DiscreteIndexedScalar(const parent_ptr& _parent, const TemporalIndexExpr indexExpr) :
    AbstractDiscreteObjectIndexed<discrete_scalar_tag>(_parent, indexExpr)
  {
  }

  virtual void accept(DiscreteExprVisitor& visitor)
  {
    visitor.visit(*this);
  }
};

class DiscreteIndexedField : public AbstractDiscreteObjectIndexed<discrete_field_tag>
{
private:
  typedef DiscreteTraits<discrete_field_tag>::expr_ptr expr_ptr;

public:
  DiscreteIndexedField(const parent_ptr& _parent, const TemporalIndexExpr indexExpr) :
    AbstractDiscreteObjectIndexed<discrete_field_tag>(_parent, indexExpr)
  {
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
};

class DiscreteIndexedOperator : public AbstractDiscreteObjectIndexed<discrete_operator_tag>
{
public:
  DiscreteIndexedOperator(const parent_ptr& _parent, const TemporalIndexExpr indexExpr) :
    AbstractDiscreteObjectIndexed<discrete_operator_tag>(_parent, indexExpr)
  {
  }

  virtual void accept(DiscreteExprVisitor& visitor)
  {
    visitor.visit(*this);
  }
};

}

}

#endif
