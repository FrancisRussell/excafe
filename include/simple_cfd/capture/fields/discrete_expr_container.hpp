#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_CONTAINER_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_CONTAINER_HPP

#include <set>
#include "discrete_traits.hpp"
#include "discrete_expr_set.hpp"
#include "temporal_index_value.hpp"
#include "function_space_expr.hpp"

namespace cfd
{

namespace detail
{

class DiscreteExprContainer
{
private:
  std::set<TemporalIndexValue*> temporalIndices;
  DiscreteExprSet<discrete_scalar_tag> scalarExpressions;
  DiscreteExprSet<discrete_field_tag> fieldExpressions;
  DiscreteExprSet<discrete_operator_tag> operatorExpressions;

public:
  bool insert(DiscreteFieldExpr& d)
  {
    return fieldExpressions.insert(d);
  }

  bool insert(ScalarExpr& d)
  {
    return scalarExpressions.insert(d);
  }

  bool insert(OperatorExpr& d)
  {
    return operatorExpressions.insert(d);
  }

  bool insert(IndexableValue<discrete_scalar_tag>& i)
  {
    return scalarExpressions.insert(i);
  }

  bool insert(IndexableValue<discrete_field_tag>& i)
  {
    return fieldExpressions.insert(i);
  }

  bool insert(IndexableValue<discrete_operator_tag>& i)
  {
    return operatorExpressions.insert(i);
  }

  bool insert(TemporalIndexValue& v)
  {
    return temporalIndices.insert(&v).second;
  }

  DiscreteExprSet<discrete_scalar_tag> getScalarExpressions() const
  {
    return scalarExpressions;
  }

  DiscreteExprSet<discrete_field_tag> getFieldExpressions() const
  {
    return fieldExpressions;
  }

  DiscreteExprSet<discrete_operator_tag> getOperatorExpressions() const
  {
    return operatorExpressions;
  }

  std::set<FunctionSpaceExpr*> getFunctionSpaces() const
  {
    std::set<FunctionSpaceExpr*> functionSpaces;

    for(DiscreteExprSet<discrete_field_tag>::expr_iter exprIter(fieldExpressions.begin_expr());
      exprIter!=fieldExpressions.end_expr(); ++exprIter)
    {
      functionSpaces.insert(&(*exprIter->getFunctionSpace()));
    }

    for(DiscreteExprSet<discrete_operator_tag>::expr_iter exprIter(operatorExpressions.begin_expr());
      exprIter!=operatorExpressions.end_expr(); ++exprIter)
    {
      functionSpaces.insert(&(*exprIter->getTrialSpace()));
      functionSpaces.insert(&(*exprIter->getTestSpace()));
    }

    return functionSpaces;
  }
};

}

}

#endif
