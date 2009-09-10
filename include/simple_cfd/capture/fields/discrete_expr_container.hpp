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
public:
  std::set<TemporalIndexValue::index_ptr> temporalIndices;
  DiscreteExprSet<discrete_scalar_tag> scalarExpressions;
  DiscreteExprSet<discrete_field_tag> fieldExpressions;
  DiscreteExprSet<discrete_operator_tag> operatorExpressions;
};

}

}

#endif
