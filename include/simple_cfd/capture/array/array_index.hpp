#ifndef SIMPLE_CFD_CAPTURE_ARRAY_ARRAY_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_ARRAY_INDEX_HPP

#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_convertible.hpp> 

namespace cfd
{

namespace detail
{

template<typename I>
class ArrayIndex
{
public:
  typedef I index_t;

private:
  BOOST_STATIC_ASSERT((boost::is_convertible<std::size_t, index_t>::value));
  std::size_t numIndices;
  std::vector<index_t> indices;

public:
  ArrayIndex(const std::size_t _numIndices) : numIndices(_numIndices), indices(numIndices)
  {
    std::fill(0, indices.begin(), indices.end());
  }

  ArrayIndex(const std::size_t _numIndices, const std::size_t* const _indices) :
    numIndices(_numIndices), indices(numIndices)
  {
    std::copy(_indices, _indices+numIndices, indices.begin());
  }

  ArrayIndex(const std::size_t _numIndices, const index_t* const _indices) :
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
