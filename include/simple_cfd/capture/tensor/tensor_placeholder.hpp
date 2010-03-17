#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_PLACEHOLDER_HPP

#include <vector>
#include "tensor_array_placeholder.hpp"
#include "scalar_placeholder.hpp"

namespace cfd
{

namespace detail
{

class TensorPlaceholder
{
private:
  TensorArrayPlaceholder tensorArray;
  ArrayIndex arrayIndex;

public:
  TensorPlaceHolder(IndexGenerator& g, const TensorArrayPlaceholder& _tensorArray) :
    tensorArray(_tensorArray), ArrayIndex(tensorArray.getArraySize())
  {
  }

  ScalarPlaceHolder operator()(const TensorIndex::constant_t i)
  {
  }
};

}

}
#endif
