#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_EXPRESSION_VALUES_SCOPE_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_EXPRESSION_VALUES_SCOPE_HPP

#include <cstddef>
#include <simple_cfd/capture/fields/discrete_traits.hpp>
#include "expression_values_map.hpp"
#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class ExpressionValuesScope
{
private:
  static const std::size_t dimension = D;
  typedef typename DiscreteValueTraits<discrete_scalar_tag, D>::value_t scalar_value_t;
  typedef typename DiscreteValueTraits<discrete_field_tag, D>::value_t field_value_t;
  typedef typename DiscreteValueTraits<discrete_operator_tag, D>::value_t operator_value_t;

  boost::shared_ptr<ExpressionValuesScope> parent;
  ExpressionValuesMap<discrete_scalar_tag, dimension> scalars;
  ExpressionValuesMap<discrete_field_tag, dimension> fields;
  ExpressionValuesMap<discrete_operator_tag, dimension> operators;

public:
  ExpressionValuesScope(const boost::shared_ptr<ExpressionValuesScope>& _parent) : parent(_parent)
  {
  }

  boost::shared_ptr<ExpressionValuesScope> getParent() const
  {
    return parent;
  }

  void completeIteration()
  {
    scalars.completeIteration();
    fields.completeIteration();
    operators.completeIteration();
  }

  scalar_value_t& getValue(ScalarExpr& e)
  {
    if (scalars.hasValue(e))
    {
      return scalars.getValue(e);
    }
    else if (parent.use_count() > 0)
    {
      return parent->getValue(e);
    }
    else
    {
      assert(false && "Unable to find value of scalar expression");
    }
  }

  field_value_t& getValue(DiscreteFieldExpr& e)
  {
    if (fields.hasValue(e))
    {
      return fields.getValue(e);
    }
    else if (parent.use_count() > 0)
    {
      return parent->getValue(e);
    }
    else
    {
      assert(false && "Unable to find value of field expression");
    }
  }

  operator_value_t& getValue(OperatorExpr& e)
  {
    if (operators.hasValue(e))
    {
      return operators.getValue(e);
    }
    else if (parent.use_count() > 0)
    {
      return parent->getValue(e);
    }
    else
    {
      assert(false && "Unable to find value of operator expression");
    }
  }

  scalar_value_t& getValue(IndexableValue<discrete_scalar_tag>& i, const signed offset)
  {
    return scalars.getValue(i, offset);
  }

  field_value_t& getValue(IndexableValue<discrete_field_tag>& i, const signed offset)
  {
    return fields.getValue(i, offset);
  }

  operator_value_t& getValue(IndexableValue<discrete_operator_tag>& i, const signed offset)
  {
    return operators.getValue(i, offset);
  }

  void setValue(ScalarExpr& e, const scalar_value_t& v)
  {
    scalars.setValue(e, v);
  }

  void setValue(DiscreteFieldExpr& e, const field_value_t& v)
  {
    fields.setValue(e, v);
  }

  void setValue(OperatorExpr& e, const operator_value_t& v)
  {
    operators.setValue(e, v);
  }

  void setValue(IndexableValue<discrete_scalar_tag>& i, const scalar_value_t& v)
  {
    scalars.setValue(i, v);
  }

  void setValue(IndexableValue<discrete_field_tag>& i, const field_value_t& v)
  {
    fields.setValue(i, v);
  }

  void setValue(IndexableValue<discrete_operator_tag>& i, const operator_value_t& v)
  {
    operators.setValue(i, v);
  }

  void calculateFinals()
  {
    assert(parent.use_count() != 0);

    const std::map<ScalarExpr*, scalar_value_t> scalarFinals = scalars.getFinals();
    parent->scalars.addMappings(scalarFinals);

    const std::map<DiscreteFieldExpr*, field_value_t> fieldFinals = fields.getFinals();
    parent->fields.addMappings(fieldFinals);

    const std::map<OperatorExpr*, operator_value_t> operatorFinals = operators.getFinals();
    parent->operators.addMappings(operatorFinals);
  }
};

}

}

#endif
