#include <set>
#include <map>
#include <cstddef>
#include <algorithm>
#include <boost/foreach.hpp>
#include <simple_cfd/capture/array/tensor_index.hpp>

namespace cfd
{

namespace detail
{

template<>
bool TensorIndex<param_tag>::isParameterised() const
{
  bool parameterised = false;

  BOOST_FOREACH(const index_t& index, indices)
  {
    parameterised |= index.isParameter();
  }

  return parameterised;
}

template<> 
std::set<TensorIndexID> TensorIndex<param_tag>::getReferencedParameters() const
{
  std::set<TensorIndexID> referenced;

  BOOST_FOREACH(const index_t& index, indices)
  {
    if (index.isParameter())
      referenced.insert(index.getParameter());
  }

  return referenced;
}

template<>
TensorIndex<param_tag> TensorIndex<param_tag>::substituteLiterals(const std::map<TensorIndexID, std::size_t>& mapping) const
{
  TensorIndex specialised(*this);

  BOOST_FOREACH(SingleIndex<TensorIndexID>& singleIndex, specialised)
  {
    if (singleIndex.isParameter())
    {
      const TensorIndexID param = singleIndex.getParameter();
      const std::map<TensorIndexID, std::size_t>::const_iterator paramIter = mapping.find(param);

      if (paramIter != mapping.end())
      {
        singleIndex = SingleIndex<TensorIndexID>(paramIter->second);
      }
    }
  }

  return specialised;
}

template<>
TensorIndex<param_tag>::TensorIndex(const std::size_t _rank, const std::size_t _dimension, 
  const TensorIndexID* const _indices) : rank(_rank), dimension(_dimension), indices(rank)
{
  std::copy(_indices, _indices+rank, indices.begin());
}

}

}
