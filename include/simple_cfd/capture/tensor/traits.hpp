#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TRAITS_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TRAITS_HPP

#include "tensor_fwd.hpp"

namespace cfd
{

namespace detail
{

template<typename object_tag>
struct IndexProperties
{
};

template<>
struct IndexProperties<array_tag>
{
  typedef ArrayIndexVariable   index_variable_t;
  typedef ArrayIndexExpression index_expression_t;
  typedef ArraySize            size_type;
};

template<>
struct IndexProperties<tensor_tag>
{
  typedef TensorIndexVariable   index_variable_t;
  typedef TensorIndexExpression index_expression_t;
  typedef TensorSize            size_type;
};

}

}

#endif
