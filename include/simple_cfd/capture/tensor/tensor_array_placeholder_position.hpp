#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_PLACEHOLDER_POSITION_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_PLACEHOLDER_POSITION_HPP

#include <cstddef>
#include "tensor_array_placeholder.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayPlaceholderPosition : public TensorArrayPlaceholder
{
public:
  TensorArrayPlaceholderPosition(const std::size_t dimension) :
    TensorArrayPlaceholder(TensorSize(1, dimension))
  {
  }
};

}

}
#endif
