#include <cstddef>
#include <cassert>
#include <map>
#include <vector>
#include <utility>
#include <boost/foreach.hpp>
#include <simple_cfd/capture/array/array_index.hpp>
#include <simple_cfd/capture/array/tensor_index.hpp>
#include <simple_cfd/capture/array/tensor_array_function_helper.hpp>
#include <simple_cfd/exception.hpp>

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

ArrayIndex<fixed_tag> TensorArrayFunctionHelper::getIndex(const std::map<ArrayIndexID, std::size_t>& parentIndices,
  const ArrayIndex<param_tag>& bindings)
{
  const ArrayIndex<param_tag> fullyBound = bindings.substituteLiterals(parentIndices);

  if (fullyBound.isParameterised())
    CFD_EXCEPTION("Attempted to create a concrete ArrayIndex without bindings for each parameter.");

  ArrayIndex<fixed_tag> result(fullyBound.numIndices());

  for(std::size_t i=0; i<fullyBound.numIndices(); ++i)
    result[i] = fullyBound[i].getConstant();

  return result;
}

TensorIndex<fixed_tag> TensorArrayFunctionHelper::getIndex(const std::map<TensorIndexID, std::size_t>& parentIndices,
  const TensorIndex<param_tag>& bindings)
{
  const TensorIndex<param_tag> fullyBound = bindings.substituteLiterals(parentIndices);

  if (fullyBound.isParameterised())
    CFD_EXCEPTION("Attempted to create a concrete ArrayIndex without bindings for each parameter.");

  TensorIndex<fixed_tag> result(fullyBound.getRank(), fullyBound.getDimension());

  for(std::size_t i=0; i<fullyBound.getRank(); ++i)
    result[i] = fullyBound[i].getConstant();

  return result;
}

 

}

}
