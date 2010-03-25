#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_PLACEHOLDER_HPP

#include <vector>
#include <cassert>
#include "tensor_fwd.hpp"
#include "index_generator.hpp"
#include "index.hpp"

namespace cfd
{

namespace detail
{

class TensorPlaceholder
{
private:
  long id;
  ArrayIndex arrayIndices;
  TensorSize tensorSize;

  void generateNewArrayIndices(IndexGenerator& g)
  {
    const ArraySize arraySize = arrayIndices.getSize();
    for(std::size_t i=0; i<arraySize.numIndices(); ++i)
    {
      arrayIndices[i] = g.newArrayIndexVariable(arraySize.getLimit(i)); 
    }
  }

public:
  TensorPlaceholder(IndexGenerator& g, const long _id, const ArraySize& _arraySize,
    const TensorSize& _tensorSize) : id(_id), arrayIndices(_arraySize), tensorSize(_tensorSize)
  {
    generateNewArrayIndices(g);
    assert(arrayIndices.allVariable());
  }

  ScalarPlaceholder operator()(const TensorIndex::constant_t i) const;
  ScalarPlaceholder operator[](const TensorIndex::constant_t i) const;
  bool operator==(const TensorPlaceholder& t) const;
  bool operator<(const TensorPlaceholder& t) const;
};

}

}
#endif
