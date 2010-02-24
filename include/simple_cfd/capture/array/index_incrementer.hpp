#ifndef SIMPLE_CFD_CAPTURE_ARRAY_INDEX_INCREMENTER_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_INDEX_INCREMENTER_HPP

#include <cstddef>
#include <map>
#include <vector>
#include <cassert>
#include <utility>
#include <boost/foreach.hpp>
#include "parameter_identifiers.hpp"

namespace cfd
{

namespace detail
{

class IndexIncrementer
{
private:
  std::vector<ArrayIndexID> arrayIndices;
  std::vector<TensorIndexID> tensorIndices;
  std::vector<std::size_t> arrayExtents;
  const std::size_t dimension;

  bool increment(std::map<ArrayIndexID, std::size_t>& arrayIndex) const
  {
    assert(arrayIndex.size() == arrayIndices.size());
    // Returns whether or not there was a wrap-around
    for(std::size_t i=arrayIndices.size()-1; i>=0; --i)
    {
      arrayIndex[arrayIndices[i]] = (arrayIndex[arrayIndices[i]] + 1) % arrayExtents[i];

      if (arrayIndex[arrayIndices[i]] != 0)
        return false;
    }
    
    return true;
  }

  bool increment(std::map<TensorIndexID, std::size_t>& tensorIndex) const
  {
    assert(tensorIndex.size() == tensorIndices.size());
    // Returns whether or not there was a wrap-around
    BOOST_REVERSE_FOREACH(const TensorIndexID& id, tensorIndices)
    {
      tensorIndex[id] = (tensorIndex[id] + 1) % dimension;

      if (tensorIndex[id] != 0)
        return false;
    }

    return true;
  }

public:
  IndexIncrementer(const std::vector<ArrayIndexID>& _arrayIndices, 
    const std::vector<TensorIndexID>& _tensorIndices,
    const std::vector<std::size_t> _arrayExtents,
    const std::size_t _dimension) :
    arrayIndices(_arrayIndices), tensorIndices(_tensorIndices), arrayExtents(_arrayExtents),
    dimension(_dimension)
  {
    assert(arrayIndices.size() == arrayExtents.size());
  }

  void zero(std::map<ArrayIndexID, std::size_t>& arrayIndex) const
  {
    arrayIndex.clear();

    BOOST_FOREACH(const ArrayIndexID& id, arrayIndices)
    {
      arrayIndex.insert(std::make_pair(id, 0));
    }
  }

  void zero(std::map<TensorIndexID, std::size_t>& tensorIndex) const
  {
    tensorIndex.clear();

    BOOST_FOREACH(const TensorIndexID& id, tensorIndices)
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

}

}

#endif
