#ifndef SIMPLE_CFD_CAPTURE_TENSOR_SCALAR_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_SCALAR_PLACEHOLDER_HPP

#include "tensor_array_placeholder.hpp"
#include "index.hpp"

namespace cfd
{

namespace detail
{

class ScalarPlaceHolder
{
private:
  TensorArrayPlaceholder tensorArray;
  ArrayIndex arrayIndex;
  TensorIndex tensorIndex;

public:
  ScalarPlaceHolder(const TensorArrayPlaceholder& _tensorArray, 
    const ArrayIndex& _arrayIndex, const TensorIndex& _tensorIndex)
  {
  }
};

}

}

#endif
