#ifndef SIMPLE_CFD_CAPTURE_DISCRETE_TRAITS_HPP
#define SIMPLE_CFD_CAPTURE_DISCRETE_TRAITS_HPP

namespace cfd
{

namespace detail
{

class discrete_scalar_tag {};
class discrete_field_tag {};
class discrete_operator_tag {};

templaye<typename T>
struct DiscreteTraits
{
};

template<>
struct DiscreteTraits<fields_scalar_tag>
{
};

template<>
struct DiscreteTraits<fields_field_tag>
{
};

template<>
struct DiscreteTraits<fields_operator_tag>
{
};

}

}

#endif
