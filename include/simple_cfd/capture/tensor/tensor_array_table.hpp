#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_HPP

#include "tensor_fwd.hpp"
#include "index_generator.hpp"
#include "array_size.hpp"
#include "tensor_size.hpp"
#include "index.hpp"

namespace cfd
{

namespace detail
{

template<typename T>
class TensorArrayTable
{
private:
  typedef T element_t;

  ArraySize arraySize;
  TensorSize tensorSize;
  ArrayIndex arrayIndices;
  TensorIndex tensorIndices;
  std::vector<element_t> table;

  std::size_t flatten(const ArrayIndex& a, const TensorIndex& t) const
  {
    return a.flatten(row_major_tag()) * tensorSize.getExtent() + t.flatten(row_major_tag());
  }

public:
  TensorArrayTable(IndexGenerator& g, const ArraySize& _arraySize, const TensorSize& _tensorSize) :
    arraySize(_arraySize), tensorSize(_tensorSize), table(arraySize.getExtent()*tensorSize.getExtent())
  {
    for(int i=0; i<arraySize.numIndices(); ++i)
    {
      arrayIndices[i] = g.newArrayIndexVariable(arraySize.getLimit(i)); 
    }

    for(int i=0; i<tensorSize.numIndices(); ++i)
    {
      tensorIndices[i] = g.newTensorIndexVariable(tensorSize.getLimit(i)); 
    }
  }

  element_t operator()(const ArrayIndex& arrayIndex, const TensorIndex& tensorIndex)
  {
    return table[flatten(arrayIndex, tensorIndex)];
  }
};

}

}

#endif
