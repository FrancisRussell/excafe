#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_HPP

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


public:
  TensorArrayTable(IndexGenerator& g, const ArraySize& _arraySize, const TensorSize& _tensorSize) :
    arraySize(_arraySize), tensorSize(_tensorSize)
  {
  }
};

}

}

#endif
