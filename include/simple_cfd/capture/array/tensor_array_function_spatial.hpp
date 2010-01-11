#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_SPATIAL_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_SPATIAL_HPP

#include <set>
#include <cstddef>
#include "array_index.hpp"
#include "parameter_identifiers.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionSpatial
{
private:
  ArrayIndex arrayExtents;
  std::size_t rank;
  std::size_t dimension;

  std::set<ArrayIndexID> virtualArrayIndices;
  std::set<TensorIndexID> virtualTensorIndices;
  
public:
  TensorArrayFunctionSpatial(const ArrayIndex& _arrayExtents, const std::size_t _rank, const std::size_t _dimension) :
    arrayExtents(_arrayExtents), rank(_rank), dimension(_dimension)
  {
  }

  TensorArrayFunctionSpatial(const ArrayIndex& _arrayExtents, const std::size_t _rank, const std::size_t _dimension,
    const std::set<ArrayIndexID>& _virtualArrayIndices, const std::set<TensorIndexID>& _virtualTensorIndices) :
    arrayExtents(_arrayExtents), rank(_rank), dimension(_dimension),
    virtualArrayIndices(_virtualArrayIndices), virtualTensorIndices(_virtualTensorIndices)
  {
  }
};

}

}

#endif
