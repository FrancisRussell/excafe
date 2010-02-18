#include <cstddef>
#include <cassert>
#include <map>
#include <vector>
#include <utility>
#include <simple_cfd/capture/array/array_index.hpp>
#include <simple_cfd/capture/array/tensor_index.hpp>
#include <simple_cfd/capture/array/tensor_array_function_helper.hpp>

namespace cfd
{

namespace detail
{

std::map<ArrayIndexID, std::size_t> TensorArrayFunctionHelper::indexToMap(const std::vector<ArrayIndexID>& indexList, 
    const ArrayIndex<fixed_tag>& index)
{
  assert(indexList.size() == index.numIndices());

  std::map<ArrayIndexID, std::size_t> mapping;

  for(std::size_t i=0; i<indexList.size(); ++i)
    mapping.insert(std::make_pair(indexList[i], index[i]));

  return mapping;
}

std::map<TensorIndexID, std::size_t> TensorArrayFunctionHelper::indexToMap(const std::vector<TensorIndexID>& indexList, 
    const TensorIndex<fixed_tag>& index)
{
  assert(indexList.size() == index.getRank());

  std::map<TensorIndexID, std::size_t> mapping;

  for(std::size_t i=0; i<indexList.size(); ++i)
    mapping.insert(std::make_pair(indexList[i], index[i]));

  return mapping;
}


}

}
