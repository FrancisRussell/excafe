#ifndef SIMPLE_CFD_CAPTURE_NUMERIC_TRAITS_HPP
#define SIMPLE_CFD_CAPTURE_NUMERIC_TRAITS_HPP

#include "numeric_fwd.hpp"

namespace cfd
{

namespace detail
{

template<typename object_tag>
struct IndexProperties
{
};

template<>
struct IndexProperties<tensor_tag>
{
  typedef TensorSize            size_type;
  typedef std::size_t           index_expression_t;
};

}

}

#endif
