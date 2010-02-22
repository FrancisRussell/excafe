#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_HELPER_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_HELPER_HPP

#include <cstddef>
#include <map>
#include <vector>
#include "array_index.hpp"
#include "tensor_index.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionHelper
{
public:
  static std::map<ArrayIndexID, std::size_t> indexToMap(const std::vector<ArrayIndexID>& indexList, 
    const ArrayIndex<fixed_tag>& index);

  static std::map<TensorIndexID, std::size_t> indexToMap(const std::vector<TensorIndexID>& indexList, 
    const TensorIndex<fixed_tag>& index);

  static ArrayIndex<fixed_tag> getIndex(const std::map<ArrayIndexID, std::size_t>& parentIndices, 
    const ArrayIndex<param_tag>& bindings);

  static TensorIndex<fixed_tag> getIndex(const std::map<TensorIndexID, std::size_t>& parentIndices, 
   const TensorIndex<param_tag>& bindings);
};

}

}

#endif
