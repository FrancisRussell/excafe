#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_HPP

#include "tensor_fwd.hpp"
#include "index_generator.hpp"
#include "array_size.hpp"
#include "tensor_size.hpp"
#include "tensor_array.hpp"
#include "index.hpp"
#include <cassert>

namespace cfd
{

namespace detail
{

template<typename T>
class TensorArrayTable : public TensorArray
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
    return ArrayIndex::flatten(a, row_major_tag()) * tensorSize.getExtent() + 
      TensorIndex::flatten(t, row_major_tag());
  }

  void generateNewArrayIndices(IndexGenerator& g)
  {
    assert(arrayIndices.getSize() == arraySize);
    for(std::size_t i=0; i<arraySize.numIndices(); ++i)
    {
      arrayIndices[i] = g.newArrayIndexVariable(arraySize.getLimit(i)); 
    }
  }

  void generateNewTensorIndices(IndexGenerator& g)
  {
    assert(tensorIndices.getSize() == tensorSize);
    for(std::size_t i=0; i<tensorSize.numIndices(); ++i)
    {
      tensorIndices[i] = g.newTensorIndexVariable(tensorSize.getLimit(i)); 
    }
  }

protected:
   ArraySize getTableArraySize() const
   {
     return arraySize;
   }

   TensorSize getTableTensorSize() const
   {
     return tensorSize;
   }

public:
  TensorArrayTable(IndexGenerator& g, const ArraySize& _arraySize, const TensorSize& _tensorSize) :
    arraySize(_arraySize), tensorSize(_tensorSize), 
    arrayIndices(arraySize), tensorIndices(tensorSize), 
    table(arraySize.getExtent()*tensorSize.getExtent())
  {
    generateNewArrayIndices(g);
    generateNewTensorIndices(g);
  }

  TensorArrayTable(IndexGenerator& g, const ArrayIndex _arrayIndices, const TensorSize& _tensorSize) :
    arraySize(_arrayIndices.getSize()), tensorSize(_tensorSize), 
    arrayIndices(_arrayIndices),
    table(arraySize.getExtent()*tensorSize.getExtent())
  {
    generateNewTensorIndices(g);
  }

  element_t& operator()(const ArrayIndex& arrayIndex, const TensorIndex& tensorIndex)
  {
    return table[flatten(arrayIndex, tensorIndex)];
  }
};

}

}

#endif
