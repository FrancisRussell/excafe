#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_PLACEHOLDER_HPP

#include <vector>
#include <cassert>
#include "tensor_array_placeholder.hpp"
#include "scalar_placeholder.hpp"
#include "index_generator.hpp"
#include "index.hpp"

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
  TensorPlaceholder(IndexGenerator& g, const TensorArrayPlaceholder& _tensorArray) :
    tensorArray(_tensorArray), arrayIndex(tensorArray.getArraySize())
  {
  }

  ScalarPlaceholder operator()(const TensorIndex::constant_t i) const
  {
    const TensorSize tensorSize = tensorArray.getTensorSize();
    const std::size_t rank = tensorSize.getRank();
    assert(rank == 1);

    TensorIndex tensorIndex(tensorSize);
    tensorIndex[0] = i;

    return ScalarPlaceholder(tensorArray, arrayIndex, tensorIndex);
  }

  ScalarPlaceholder operator[](const TensorIndex::constant_t i) const
  {
    return (*this)(i);
  }
};

}

}
#endif
