#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_EXPRESSION_VALUES_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_EXPRESSION_VALUES_HPP

#include <cstddef>
#include <simple_cfd/capture/fields/discrete_traits.hpp>
#include "expression_values_typed.hpp"

namespace cfd
{

namespace detail
{

template<std::size_t D>
class ExpressionValues
{
private:
  static const std::size_t dimension = D;

  ExpressionValuesTyped<discrete_scalar_tag, dimension> scalars;
  ExpressionValuesTyped<discrete_field_tag, dimension> fields;
  ExpressionValuesTyped<discrete_operator_tag, dimension> operators;

};

}

}

#endif
