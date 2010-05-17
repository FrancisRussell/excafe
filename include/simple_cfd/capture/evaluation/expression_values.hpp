#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_EXPRESSION_VALUES_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_EXPRESSION_VALUES_HPP

#include <cstddef>
#include <cassert>
#include <boost/shared_ptr.hpp>
#include <simple_cfd/capture/fields/discrete_traits.hpp>
#include "expression_values_scope.hpp"

namespace cfd
{

namespace detail
{

template<std::size_t D>
class ExpressionValues
{
private:
  static const std::size_t dimension = D;
  typedef typename DiscreteValueTraits<discrete_scalar_tag, D>::value_t scalar_value_t;
  typedef typename DiscreteValueTraits<discrete_field_tag, D>::value_t field_value_t;
  typedef typename DiscreteValueTraits<discrete_operator_tag, D>::value_t operator_value_t;

  boost::shared_ptr< ExpressionValuesScope<dimension> > values;

public:
  ExpressionValues()
  {
  }

  void enterScope()
  {
    values = boost::shared_ptr< ExpressionValuesScope<dimension> >(new ExpressionValuesScope<dimension>(values));
  }

  void exitScope()
  {
    assert(values.use_count() != 0);
    values = values->getParent();
  }

  void calculateInitials(const TemporalIndexValue& loopIndex)
  {
    values->calculateInitials(loopIndex);
  }

  void calculateFinals()
  {
    values->calculateFinals();
  }

  void completeIteration()
  {
    assert(!values.use_count() == 0);
    values->completeIteration();
  }

  bool isGlobalScope() const
  {
    return values->getParent().use_count() == 0;
  }

  scalar_value_t& getValue(ScalarExpr& e)
  {
    return values->getValue(e);
  }

  field_value_t& getValue(DiscreteFieldExpr& e)
  {
    return values->getValue(e);
  }

  operator_value_t& getValue(OperatorExpr& e)
  {
    return values->getValue(e);
  }

  bool hasValue(ScalarExpr& e) const
  {
    return values->hasValue(e);
  }

  bool hasValue(DiscreteFieldExpr& e) const
  {
    return values->hasValue(e);
  }

  bool getValue(OperatorExpr& e) const
  {
    return values->hasValue(e);
  }

  bool hasValue(IndexableValue<discrete_scalar_tag>& i, const signed offset) const
  {
    return values->hasValue(i, offset);
  }

  bool hasValue(IndexableValue<discrete_field_tag>& i, const signed offset) const
  {
    return values->hasValue(i, offset);
  }

  bool hasValue(IndexableValue<discrete_operator_tag>& i, const signed offset) const
  {
    return values->hasValue(i, offset);
  }

  scalar_value_t& getValue(IndexableValue<discrete_scalar_tag>& i, const signed offset)
  {
    return values->getValue(i, offset);
  }

  field_value_t& getValue(IndexableValue<discrete_field_tag>& i, const int offset)
  {
    return values->getValue(i, offset);
  }

  operator_value_t& getValue(IndexableValue<discrete_operator_tag>& i, const int offset)
  {
    return values->getValue(i, offset);
  }

  void setValue(ScalarExpr& e, const scalar_value_t& v)
  {
    values->setValue(e, v);
  }

  void setValue(DiscreteFieldExpr& e, const field_value_t& v)
  {
    values->setValue(e, v);
  }

  void setValue(OperatorExpr& e, const operator_value_t& v)
  {
    values->setValue(e, v);
  }

  void  setValue(IndexableValue<discrete_scalar_tag>& i, const scalar_value_t& v, const signed offset)
  {
    values->setValue(i, v, offset);
  }

  void  setValue(IndexableValue<discrete_field_tag>& i, const field_value_t& v, const signed offset)
  {
    values->setValue(i, v, offset);
  }

  void  setValue(IndexableValue<discrete_operator_tag>& i, const operator_value_t& v, const signed offset)
  {
    values->setValue(i, v, offset);
  }
};

}

}

#endif
