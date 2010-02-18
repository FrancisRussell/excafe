#include <set>
#include <map>
#include <cstddef>
#include <algorithm>
#include <boost/foreach.hpp>
#include <simple_cfd/capture/array/array_index.hpp>

namespace cfd
{

namespace detail
{

template<>
bool ArrayIndex<param_tag>::isParameterised() const
{
  bool parameterised = false;

  BOOST_FOREACH(const index_t& index, indices)
  {
    parameterised |= index.isParameter();
  }
    
  return parameterised;
}

template<> 
std::set<ArrayIndexID> ArrayIndex<param_tag>::getReferencedParameters() const
{
  std::set<ArrayIndexID> referenced;

  BOOST_FOREACH(const index_t& index, indices)
  {
    if (index.isParameter())
      referenced.insert(index.getParameter());
  }

  return referenced;
}

template<>
ArrayIndex<param_tag> ArrayIndex<param_tag>::substituteLiterals(const std::map<ArrayIndexID, std::size_t>& mapping) const
{
  ArrayIndex specialised(*this);

  BOOST_FOREACH(SingleIndex<ArrayIndexID>& singleIndex, specialised)
  {
    if (singleIndex.isParameter())
    {
      const ArrayIndexID param = singleIndex.getParameter();
      const std::map<ArrayIndexID, std::size_t>::const_iterator paramIter = mapping.find(param);

      if (paramIter != mapping.end())
      {
        singleIndex = SingleIndex<ArrayIndexID>(paramIter->second);
      }
    }
  }

  return specialised;
}

template<>
ArrayIndex<param_tag>::ArrayIndex(const std::size_t _numIndices, const ArrayIndexID* const _indices) : indices(_numIndices)
{
  std::copy(_indices, _indices+_numIndices, indices.begin());
}

}

}
