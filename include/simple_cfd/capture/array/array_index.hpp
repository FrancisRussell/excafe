#ifndef SIMPLE_CFD_CAPTURE_ARRAY_ARRAY_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_ARRAY_INDEX_HPP

#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>

namespace cfd
{

namespace detail
{

class ArrayIndex
{
private:
  std::size_t numIndices;
  std::vector<std::size_t> indices;

public:
  ArrayIndex(const std::size_t _numIndices) : numIndices(_numIndices)
  {
  }

  ArrayIndex(const std::size_t _numIndices, const std::size_t* const _indices) :
    numIndices(_numIndices), indices(numIndices)
  {
    std::copy(_indices, _indices+numIndices, indices.begin());
  }

  std::size_t& operator[](const std::size_t index)
  {
    assert(index < indices.size());
    return indices[index];
  }

  std::size_t operator[](const std::size_t index) const
  {
    assert(index < indices.size());
    return indices[index];
  }

  std::size_t getNumIndices() const
  {
    return numIndices;
  }

  std::size_t getDimension(const std::size_t i) const
  {
    assert(i<indices.size());
    return indices[i];
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
