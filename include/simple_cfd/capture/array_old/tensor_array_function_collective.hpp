#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_COLLECTIVE_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_COLLECTIVE_HPP

#include <cstddef>
#include <vector>
#include <map>
#include <utility>
#include <boost/tuple/tuple.hpp>
#include "tensor_function.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "tensor_array_function_call.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionCollective : public TensorFunction
{
protected:
  typedef TensorArrayFunctionCall call_t;
  const ArrayIndex<fixed_tag> arrayExtents;
  const std::size_t rank;
  const std::size_t dimension;

  std::vector<ArrayIndexID> arrayIndexParameters;
  std::vector<TensorIndexID> tensorIndexParameters;

  std::vector<ArrayIndexID> unseenArrayIndexParameters;
  std::vector<TensorIndexID> unseenTensorIndexParameters;
  ArrayIndex<fixed_tag> unseenArrayExtents;

  std::vector<call_t> operands;

public:
  TensorArrayFunctionCollective(const ArrayIndex<fixed_tag>& _arrayExtents, const std::size_t _rank, 
    const std::size_t _dimension) : arrayExtents(_arrayExtents), rank(_rank), dimension(_dimension)
  {
    for(std::size_t i=0; i<arrayExtents.numIndices(); ++i)
       arrayIndexParameters.push_back(ArrayIndexID(i));

    for(std::size_t i=0; i<rank; ++i)
       tensorIndexParameters.push_back(TensorIndexID(i));
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
    return arrayIndexParameters.size();
  }

  virtual std::size_t getArrayDimension(const std::size_t index) const
  {
    return arrayExtents[index];
  }

  virtual ArrayIndex<fixed_tag> getArrayExtent() const
  {
    return arrayExtents;
  }

  TensorIndexID newTensorIndex()
  {
    const TensorIndexID index(tensorIndexParameters.size() + unseenTensorIndexParameters.size());
    unseenTensorIndexParameters.push_back(index);
    return index;
  }

  std::vector<TensorIndexID> newTensorIndices(const std::size_t count)
  {
    std::vector<TensorIndexID> paramVector;

    for(std::size_t i=0; i<count; ++i);
      paramVector.push_back(newTensorIndex());

    return paramVector;
  }

  ArrayIndexID newArrayIndex(const std::size_t extent)
  {
    const ArrayIndexID index(arrayIndexParameters.size() + unseenArrayIndexParameters.size());
    unseenArrayIndexParameters.push_back(index);
    unseenArrayExtents.append(extent);
    return index;
  }

  void addTerm(const ArrayIndex<param_tag>& arrayIndex, const TensorIndex<param_tag>& tensorIndex,
    const TensorFunction::ref function)
  {
    operands.push_back(call_t(arrayIndex, tensorIndex, function));
  }

  ArrayIndex<param_tag> getIdentityArrayIndex() const
  {
    return ArrayIndex<param_tag>(arrayIndexParameters.size(), &arrayIndexParameters[0]);
  }

  TensorIndex<param_tag> getIdentityTensorIndex() const
  {
    return TensorIndex<param_tag>(tensorIndexParameters.size(), dimension, &tensorIndexParameters[0]);
  }
};

}

}

#endif
