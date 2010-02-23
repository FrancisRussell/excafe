#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_HPP

#include <cstddef>
#include <numeric>
#include <functional>
#include <algorithm>
#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <iterator>
#include <simple_cfd/exception.hpp>
#include <boost/foreach.hpp>
#include "tensor_function.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "array_traits.hpp"
#include "tensor_array_function_helper.hpp"

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

  class IndexIncrementer
  {
  private:
    std::vector< std::pair<ArrayIndexID, std::size_t> > realArrayIndexExtents;
    std::vector<TensorIndexID> realTensorIndices;
    const std::size_t dimension;

    bool increment(std::map<ArrayIndexID, std::size_t>& arrayIndex) const
    {
      assert(arrayIndex.size() == realArrayIndexExtents.size());
      // Returns whether or not there was a wrap-around
      typedef std::pair<ArrayIndexID, std::size_t> extent_t;
      BOOST_REVERSE_FOREACH(const extent_t& extent, realArrayIndexExtents)
      {
        arrayIndex[extent.first] = (arrayIndex[extent.first] + 1) % extent.second;

        if (arrayIndex[extent.first] != 0)
          return false;
      }
      
      return true;
    }

    bool increment(std::map<TensorIndexID, std::size_t>& tensorIndex) const
    {
      assert(tensorIndex.size() == realTensorIndices.size());
      // Returns whether or not there was a wrap-around
      BOOST_REVERSE_FOREACH(const TensorIndexID& id, realTensorIndices)
      {
        tensorIndex[id] = (tensorIndex[id] + 1) % dimension;

        if (tensorIndex[id] != 0)
          return false;
      }

      return true;
    }

  public:
    IndexIncrementer(const TensorArrayFunction& parent) : dimension(parent.dimension)
    {
      for(std::size_t i=0; i<parent.arrayExtents.numIndices(); ++i)
      {
        if (!parent.isVirtual(parent.arrayIndexParameters[i]))
          realArrayIndexExtents.push_back(std::make_pair(parent.arrayIndexParameters[i], parent.arrayExtents[i]));
      }

      for(std::size_t i=0; i<parent.rank; ++i)
      {
        if (!parent.isVirtual(parent.tensorIndexParameters[i]))
          realTensorIndices.push_back(parent.tensorIndexParameters[i]);
      }
    }

    void zero(std::map<ArrayIndexID, std::size_t>& arrayIndex) const
    {
      arrayIndex.clear();

      typedef std::pair<ArrayIndexID, std::size_t> extent_t;
      BOOST_FOREACH(const extent_t& extent, realArrayIndexExtents)
      {
        arrayIndex.insert(std::make_pair(extent.first, 0));
      }
    }

    void zero(std::map<TensorIndexID, std::size_t>& tensorIndex) const
    {
      tensorIndex.clear();

      BOOST_FOREACH(const TensorIndexID& id, realTensorIndices)
      {
        tensorIndex.insert(std::make_pair(id, 0));
      }
    }

    bool increment(std::map<ArrayIndexID, std::size_t>& arrayIndex, 
      std::map<TensorIndexID, std::size_t>& tensorIndex) const
    {
      if (increment(tensorIndex))
        return increment(arrayIndex);
      else
        return false;
    }
  };

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
