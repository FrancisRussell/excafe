#ifndef SIMPLE_CFD_CAPTURE_ARRAY_ARRAY_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_ARRAY_INDEX_HPP

#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_convertible.hpp> 
#include "single_index.hpp"

namespace cfd
{

namespace detail
{

class ArrayIndex
{
public:
  typedef ArraySingleIndex index_t;

private:
  std::size_t numIndices;
  std::vector<index_t> indices;

public:
  ArrayIndex(const std::size_t _numIndices) : numIndices(_numIndices), indices(numIndices)
  {
    std::fill(indices.begin(), indices.end(), 0);
  }

  ArrayIndex(const std::size_t _numIndices, const index_t::constant_t* const _indices) :
    numIndices(_numIndices), indices(numIndices)
  {
    std::copy(_indices, _indices+numIndices, indices.begin());
  }

  ArrayIndex(const std::size_t _numIndices, const index_t::parameter_t* const _indices) :
    numIndices(_numIndices), indices(numIndices)
  {
    std::copy(_indices, _indices+numIndices, indices.begin());
  }

  index_t& operator[](const std::size_t index)
  {
    assert(index < indices.size());
    return indices[index];
  }

  index_t operator[](const std::size_t index) const
  {
    assert(index < indices.size());
    return indices[index];
  }

  std::size_t getNumIndices() const
  {
    return numIndices;
  }

  bool operator==(const ArrayIndex& i) const
  {
    return numIndices == i.numIndices &&
           indices == i.indices;
  }

  bool operator<(const ArrayIndex& i) const
  {
    if (numIndices < i.numIndices) return true;
    if (numIndices  == i.numIndices && indices < i.indices) return true;
    return false;
  }

  void swap(ArrayIndex& i)
  {
    std::swap(numIndices, i.numIndices);
    std::swap(indices, i.indices);
  }
};

}

}

#endif
