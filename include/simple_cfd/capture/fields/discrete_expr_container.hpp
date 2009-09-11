#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_CONTAINER_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_CONTAINER_HPP

#include <set>
#include "discrete_traits.hpp"
#include "discrete_expr_set.hpp"
#include "temporal_index_value.hpp"

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
};

}

}

#endif
