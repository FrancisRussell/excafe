#ifndef EXCAFE_DISCRETE_VALUE_TRAITS_HPP
#define EXCAFE_DISCRETE_VALUE_TRAITS_HPP

#include <cstddef>
#include <excafe/capture/fields/discrete_traits.hpp>
#include "discrete_field.hpp"
#include "discrete_operator.hpp"

namespace excafe
{

namespace detail
{

// Map value tags to value classes

template<typename T, std::size_t D>
struct DiscreteValueTraits
{
};

template<std::size_t D>
struct DiscreteValueTraits<discrete_scalar_tag, D>
{
  typedef double value_t;
};

template<std::size_t D>
struct DiscreteValueTraits<discrete_field_tag, D>
{
  typedef DiscreteField<D> value_t;
};

template<std::size_t D>
struct DiscreteValueTraits<discrete_operator_tag, D>
{
  typedef DiscreteOperator<D> value_t;
};

}

}

#endif
