#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_HPP

#include <cstddef>
#include <functional>
#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <simple_cfd/exception.hpp>
#include "tensor_function.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "array_traits.hpp"
#include "tensor_array_function_helper.hpp"
#include "index_incrementer.hpp"

namespace cfd
{

namespace detail
{

template<typename E>
class TensorArrayFunction : public TensorFunction
{
protected:
  typedef E element_t;

  const ArrayIndex<fixed_tag> arrayExtents;
  const std::size_t rank;
  const std::size_t dimension;

  std::vector<ArrayIndexID> arrayIndexParameters;
  std::vector<TensorIndexID> tensorIndexParameters;
  std::set<ArrayIndexID> arrayVirtualParameters;
  std::set<TensorIndexID> tensorVirtualParameters;
  std::vector<element_t> values;

  std::size_t getInternalRank() const
  {
    return rank - tensorVirtualParameters.size();
  }

  std::size_t getNumInternalArrayIndices() const
  {
    return arrayExtents.numIndices() - arrayVirtualParameters.size();
  }

  std::size_t tensorExtent(const bool real) const
  {
    std::size_t result = 1;

    for(std::size_t i=0; i<rank; ++i)
    {
      if (!real || tensorVirtualParameters.find(tensorIndexParameters[i]) == tensorVirtualParameters.end())
        result *= dimension;
    }

    return result;
  }

  bool isVirtual(const ArrayIndexID& index) const
  {
    return arrayVirtualParameters.find(index) != arrayVirtualParameters.end();
  }

  bool isVirtual(const TensorIndexID& index) const
  {
    return tensorVirtualParameters.find(index) != tensorVirtualParameters.end();
  }

  std::size_t tensorExtentVirtual() const
  {
    return tensorExtent(false);
  }

  std::size_t tensorExtentReal() const
  {
    return tensorExtent(true);
  }

  std::size_t arrayExtent(const bool real) const
  {
    std::size_t result = 1;

    for(std::size_t i=0; i<arrayExtents.numIndices(); ++i)
    {
      if (!real || arrayVirtualParameters.find(arrayIndexParameters[i]) == arrayVirtualParameters.end())
        result *= arrayExtents[i];
    }

    return result;
  }

  std::size_t arrayExtentVirtual() const
  {
    return arrayExtent(false);
  }

  std::size_t arrayExtentReal() const
  {
    return arrayExtent(true);
  }

  std::size_t extentVirtual() const
  {
    return arrayExtentVirtual() * tensorExtentVirtual();
  }

  std::size_t extentReal() const
  {
    return arrayExtentReal() * tensorExtentReal();
  }

  //TODO: too much code redundancy here

  std::size_t flattenReal(const ArrayIndex<fixed_tag>& arrayIndex) const
  {
    assert(arrayIndex.numIndices() == getNumInternalArrayIndices());
    
    std::size_t offset = 0;
    std::size_t multiplier = 1;
    std::size_t realIndex = arrayIndex.numIndices() - 1;

    for(int i=arrayExtents.numIndices()-1; i>=0; --i)
    {
      if (arrayVirtualParameters.find(arrayIndexParameters[i]) == arrayVirtualParameters.end())
      {
        offset += arrayIndex[realIndex] * multiplier;
        multiplier *= arrayExtents[i];
        --realIndex;
      }
    }

    assert(realIndex = 0);
    return offset;
  }

  std::size_t flattenReal(const TensorIndex<fixed_tag>& tensorIndex) const
  {
    assert(tensorIndex.getRank() == getInternalRank());
    assert(tensorIndex.getDimension() == dimension);

    std::size_t offset = 0;
    std::size_t multiplier = 1;

    for(int i=getInternalRank()-1; i>=0; --i)
    {
      if (tensorVirtualParameters.find(tensorIndexParameters[i]) == tensorVirtualParameters.end())
      {
        offset += tensorIndex[i] * multiplier;
        multiplier *= dimension;
      }
    }

    return offset;
  }

  std::size_t flattenReal(const std::map<ArrayIndexID, std::size_t>& arrayIndex) const
  {
    assert(arrayIndex.size() >= getNumInternalArrayIndices());
    
    std::size_t offset = 0;
    std::size_t multiplier = 1;

    for(int i=arrayExtents.numIndices()-1; i>=0; --i)
    {
      if (arrayVirtualParameters.find(arrayIndexParameters[i]) == arrayVirtualParameters.end())
      {
        const std::map<ArrayIndexID, std::size_t>::const_iterator arrIndexIter = arrayIndex.find(arrayIndexParameters[i]);
        if (arrIndexIter == arrayIndex.end())
        {
          CFD_EXCEPTION("Missing array index value when flattening.");
        }
        else
        {
          offset += arrIndexIter->second * multiplier;
          multiplier *= arrayExtents[i];
        }
      }
    }

    return offset;
  }

  std::size_t flattenReal(const std::map<TensorIndexID, std::size_t>& tensorIndex) const
  {
    assert(tensorIndex.size() >= getInternalRank());
    
    std::size_t offset = 0;
    std::size_t multiplier = 1;

    for(int i=rank-1; i>=0; --i)
    {
      if (tensorVirtualParameters.find(tensorIndexParameters[i]) == tensorVirtualParameters.end())
      {
        const std::map<TensorIndexID, std::size_t>::const_iterator tensorIndexIter = tensorIndex.find(tensorIndexParameters[i]);
        if (tensorIndexIter == tensorIndex.end())
        {
          CFD_EXCEPTION("Missing tensor index value when flattening.");
        }
        else
        {
          offset += tensorIndexIter->second * multiplier;
          multiplier *= dimension;
        }
      }
    }

    return offset;
  }

  std::size_t flattenReal(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex) const
  {
    return flattenReal(arrayIndex) * tensorExtentReal() + flattenReal(tensorIndex);
  }

  std::size_t flattenReal(const std::map<ArrayIndexID, std::size_t>& arrayIndex, 
    const std::map<TensorIndexID, std::size_t>& tensorIndex) const
  {
    return flattenReal(arrayIndex) * tensorExtentReal() + flattenReal(tensorIndex);
  }

  element_t& operator()(const std::map<ArrayIndexID, std::size_t>& arrayIndex, 
    const std::map<TensorIndexID, std::size_t>& tensorIndex)
  {
    const std::size_t offset = flattenReal(arrayIndex, tensorIndex);
    return values[offset];
  }

  const element_t operator()(const std::map<ArrayIndexID, std::size_t>& arrayIndex, 
    const std::map<TensorIndexID, std::size_t>& tensorIndex) const
  {
    const std::size_t offset = flattenReal(arrayIndex, tensorIndex);
    return values[offset];
  }

  void accept(TensorArrayFunctionVisitor<element_t>& visitor)
  {
    std::vector<ArrayIndexID> realArrayIndices;
    std::vector<std::size_t> realArrayExtents;

    for(std::size_t i=0; i<arrayExtents.numIndices(); ++i)
    {
      if (!isVirtual(arrayIndexParameters[i]))
      {
        realArrayIndices.push_back(arrayIndexParameters[i]);
        realArrayExtents.push_back(arrayExtents[i]);
      }
    }

    std::vector<TensorIndexID> realTensorIndices;

    for(std::size_t i=0; i<rank; ++i)
    {
      if (!isVirtual(tensorIndexParameters[i]))
        realTensorIndices.push_back(tensorIndexParameters[i]);
    }

    IndexIncrementer incrementer(realArrayIndices, realTensorIndices, realArrayExtents, dimension);
    std::map<ArrayIndexID, std::size_t> arrayIndex; 
    std::map<TensorIndexID, std::size_t> tensorIndex;

    incrementer.zero(arrayIndex);
    incrementer.zero(tensorIndex);

    do
    {
      visitor.visit(*this, arrayIndex, tensorIndex, (*this)(arrayIndex, tensorIndex));
    }
    while(!incrementer.increment(arrayIndex, tensorIndex));
  }

  TensorArrayFunction(const ArrayIndex<fixed_tag>& _arrayExtents, const std::size_t _rank, 
    const std::size_t _dimension) :
    arrayExtents(_arrayExtents), rank(_rank), dimension(_dimension), values(extentReal())
  {
    for(std::size_t i=0; i<arrayExtents.numIndices(); ++i)
       arrayIndexParameters.push_back(ArrayIndexID(i));

    for(std::size_t i=0; i<rank; ++i)
       tensorIndexParameters.push_back(TensorIndexID(i));
  }

  TensorArrayFunction(const ArrayIndex<fixed_tag>& _arrayExtents, const std::size_t _rank, 
    const std::size_t _dimension, const std::vector<ArrayIndexID>& _arrayIndexParameters,
    const std::vector<TensorIndexID>& _tensorIndexParameters) :
    arrayExtents(_arrayExtents), rank(_rank), dimension(_dimension), 
    arrayIndexParameters(_arrayIndexParameters), tensorIndexParameters(_tensorIndexParameters),
    values(extentReal())
  {
  }

public:
  virtual std::size_t getTensorRank() const
  {
    return rank;
  }

  virtual std::size_t getTensorDimension() const
  {
    return dimension;
  }

  virtual std::size_t numArrayIndices() const
  {
    return arrayExtents.numIndices();
  }

  virtual std::size_t getArrayDimension(const std::size_t index) const
  {
    return arrayExtents[index];
  }

  virtual ArrayIndex<fixed_tag> getArrayExtent() const
  {
    return arrayExtents;
  }

  ArrayIndex<param_tag> getIdentityArrayIndex() const
  {
    assert(arrayExtents.numIndices() == arrayIndexParameters.size());
    return ArrayIndex<param_tag>(arrayIndexParameters.size(), &arrayIndexParameters[0]);
  }

  TensorIndex<param_tag> getIdentityTensorIndex() const
  {
    assert(rank == tensorIndexParameters.size());
    return TensorIndex<param_tag>(rank, dimension, &tensorIndexParameters[0]);
  }

  element_t& operator()(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex)
  {
    const std::size_t offset = flattenReal(arrayIndex, tensorIndex);
    return values[offset];
  }

  const element_t operator()(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex) const
  {
    const std::size_t offset = flattenReal(arrayIndex, tensorIndex);
    return values[offset];
  }
};

}

}

#endif
