#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_HPP

#include "tensor_fwd.hpp"
#include "index_generator.hpp"
#include "array_size.hpp"
#include "tensor_size.hpp"
#include "tensor_array.hpp"
#include "index.hpp"
#include <cassert>
#include <vector>

namespace cfd
{

namespace detail
{

template<typename T>
class TensorArrayTable : public TensorArray
{
private:
  typedef T element_t;

  ArrayIndex tableArrayIndices;
  TensorIndex tableTensorIndices;
  TensorIndex additionalTensorIndices;
  std::vector<element_t> table;

  std::size_t flatten(const ArrayIndex& a, const TensorIndex& t) const
  {
    return ArrayIndex::flatten(a, row_major_tag()) * tableTensorSize().getExtent() + 
      TensorIndex::flatten(t, row_major_tag());
  }

  void generateNewArrayIndices(IndexGenerator& g)
  {
    for(std::size_t i=0; i<tableArrayIndices.numIndices(); ++i)
    {
      tableArrayIndices[i] = g.newArrayIndexVariable(tableArraySize().getLimit(i)); 
    }
  }

  void generateNewTensorIndices(IndexGenerator& g)
  {
    for(std::size_t i=0; i<tableTensorIndices.numIndices(); ++i)
    {
      tableTensorIndices[i] = g.newTensorIndexVariable(tableTensorSize().getLimit(i)); 
    }
  }

  std::size_t tableSize() const
  {
    return tableArraySize().getExtent() * tableTensorSize().getExtent(); 
  }

protected:
  ArraySize tableArraySize() const
  {
    return tableArrayIndices.getSize();
  }

  TensorSize tableTensorSize() const
  {
    return tableTensorIndices.getSize();
  }

public:
   typedef typename std::vector<element_t>::iterator iterator;
   typedef typename std::vector<element_t>::const_iterator const_iterator;

  TensorArrayTable(IndexGenerator& g, const ArraySize& _arraySize, const TensorSize& _tensorSize) :
    tableArrayIndices(_arraySize), tableTensorIndices(_tensorSize), 
    additionalTensorIndices(TensorSize(0, _tensorSize.getDimension())),
    table(tableSize())
  {
    generateNewArrayIndices(g);
    generateNewTensorIndices(g);
  }

  TensorArrayTable(IndexGenerator& g, const ArrayIndex _arrayIndices, const TensorSize& _tensorSize) :
    tableArrayIndices(_arrayIndices), tableTensorIndices(_tensorSize),
    additionalTensorIndices(TensorSize(0, _tensorSize.getDimension())),
    table(tableSize())
  {
    generateNewTensorIndices(g);
  }

  virtual TensorSize getTensorSize() const
  {
    return TensorSize(tableTensorSize().getRank() + additionalTensorIndices.getSize().getRank(), tableTensorSize().getDimension()); 
  }  

  void appendAdditionalTensorIndices(IndexGenerator& g, const std::size_t count)
  {
    for(std::size_t i=0; i<count; ++i)
      additionalTensorIndices = additionalTensorIndices.append(g.newTensorIndexVariable(tableTensorSize().getDimension()));
  }

  TensorIndex getAdditionalIndices() const
  {
    return additionalTensorIndices;
  }

  iterator begin()
  {
    return table.begin();
  }

  iterator end()
  {
    return table.end();
  }

  const_iterator begin() const
  {
    return table.begin();
  }

  const_iterator end() const
  {
    return table.end();
  }

  element_t& operator()(const ArrayIndex& arrayIndex, const TensorIndex& tensorIndex)
  {
    return table[flatten(arrayIndex, tensorIndex)];
  }
};

}

}

#endif
