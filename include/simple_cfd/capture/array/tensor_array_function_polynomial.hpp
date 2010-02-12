#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_POLYNOMIAL_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_POLYNOMIAL_HPP

#include <cstddef>
#include <numeric>
#include <functional>
#include <vector>
#include <simple_cfd/exception.hpp>
#include "tensor_function.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "array_traits.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionPolynomial : public TensorFunction
{
private:
  typedef TensorFunction::polynomial_t polynomial_t;

  const ArrayIndex<fixed_tag> arrayExtents;
  const std::size_t rank;
  const std::size_t dimension;

  std::vector<ArrayIndexID> arrayIndexParameters;
  std::vector<TensorIndexID> tensorIndexParameters;
  std::set<ArrayIndexID> arrayVirtualParameters;
  std::set<TensorIndexID> tensorVirtualParameters;
  std::vector<polynomial_t> values;

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
      if (!real || tensorVirtualParameters.find(tensorIndexParameters[i]) != tensorVirtualParameters.end())
        result *= dimension;
    }

    return result;
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
      if (!real || arrayVirtualParameters.find(arrayIndexParameters[i]) != arrayVirtualParameters.end())
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

  std::size_t flattenVirtual(const ArrayIndex<fixed_tag>& arrayIndex) const
  {
    assert(arrayIndex.numIndices() ==  arrayExtents.numIndices());
    
    std::size_t offset = 0;
    std::size_t multiplier = 1;

    for(int i=arrayIndex.numIndices()-1; i>=0; --i)
    {
      offset += arrayIndex[i] * multiplier;
      multiplier *= arrayExtents[i];
    }

    return offset;
  }

  std::size_t flattenVirtual(const TensorIndex<fixed_tag>& tensorIndex) const
  {
    assert(tensorIndex.getRank() == rank);
    assert(tensorIndex.getDimension() == dimension);

    std::size_t offset = 0;
    std::size_t multiplier = 1;

    for(int i=rank-1; i>=0; --i)
    {
      offset += tensorIndex[i] * multiplier;
      multiplier *= dimension;
    }

    return offset;
  }

  std::size_t flattenVirtual(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex) const
  {
    return flattenVirtual(arrayIndex) * tensorExtentVirtual() + flattenVirtual(tensorIndex);
  }

  std::size_t flattenReal(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex) const
  {
    //FIXME: implement me!
  }

  TensorFunction::ref expand(const std::set<ArrayIndexID>& arrayIndicesToExpand, const
    std::set<TensorIndexID>& tensorIndicesToExpand) const
  {
    //FIXME: implement me!
  }

  TensorArrayFunctionPolynomial(const ArrayIndex<fixed_tag>& _arrayExtents, 
    const std::size_t _rank, const std::size_t _dimension, 
    const std::vector<ArrayIndexID>& _arrayIndices, const std::vector<TensorIndexID>& _tensorIndices,
    const std::set<ArrayIndexID>& _virtualArrayIndices, const std::set<TensorIndexID>& _virtualTensorIndices) :
    arrayExtents(_arrayExtents), rank(_rank), dimension(_dimension), 
    arrayIndexParameters(_arrayIndices), tensorIndexParameters(_tensorIndices), 
    arrayVirtualParameters(_virtualArrayIndices), tensorVirtualParameters(_virtualTensorIndices), 
    values(extentReal())
    {
    }

public:
  TensorArrayFunctionPolynomial(const ArrayIndex<fixed_tag>& _arrayExtents, const std::size_t _rank, 
    const std::size_t _dimension) :
    arrayExtents(_arrayExtents), rank(_rank), dimension(_dimension), values(extentReal())
  {
    // FIXME: create indices
  }

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

  virtual TensorFunction::ref differentiate(const ScalarReference& reference) const
  {
    if (reference.isBound()) 
        CFD_EXCEPTION("Cannot differentiate with respect to bound reference.");

    if (reference.isParameterised()) 
      CFD_EXCEPTION("Cannot differentiate with respect to parameterised reference.");

    // FIXME: implement me!
    // 1. Check for polynomials with unbound references, that are parameterised
    // 2. If the parameter is known, substitute it, else expand out the parameter
    // 3. Conventional differentiation step on the entire tensor function
  }

  polynomial_t& operator()(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex)
  {
    const std::size_t offset = flattenReal(arrayIndex, tensorIndex);
    return values[offset];
  }
};

}

}

#endif
