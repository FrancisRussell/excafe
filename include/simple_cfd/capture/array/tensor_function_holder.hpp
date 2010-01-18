#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_FUNCTION_HOLDER_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_FUNCTION_HOLDER_HPP

#include "tensor_function.hpp"

namespace cfd
{

namespace detail
{

class TensorFunctionHolder
{
private:
  typedef TensorFunction::function_ptr function_ptr;
  function_ptr function;

public:
  std::size_t getTensorRank() const
  {
    return function->getTensorRank();
  }

  std::size_t getTensorDimension() const
  {
    return function->getTensorDimension();
  }

  std::size_t numArrayIndices() const
  {
    return function->numArrayIndices();
  }

  std::size_t getArrayDimension(const std::size_t index) const
  {
    return function->getArrayDimension(index);
  }
};

}

}

#endif
