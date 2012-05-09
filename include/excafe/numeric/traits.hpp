#ifndef EXCAFE_CAPTURE_NUMERIC_TRAITS_HPP
#define EXCAFE_CAPTURE_NUMERIC_TRAITS_HPP

#include "numeric_fwd.hpp"

namespace excafe
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
