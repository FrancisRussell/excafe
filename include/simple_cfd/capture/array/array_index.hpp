#ifndef SIMPLE_CFD_CAPTURE_ARRAY_ARRAY_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_ARRAY_INDEX_HPP

#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>
#include <boost/mpl/map.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/at.hpp>
#include "single_index.hpp"
#include "parameter_identifiers.hpp"
#include "array_traits.hpp"

namespace cfd
{

namespace detail
{

template<typename I>
class ArrayIndex
{
public:
  typedef std::size_t constant_t;
  typedef SingleIndex<ArrayIndexID> parameter_t;

private:
  typedef boost::mpl::map<
      boost::mpl::pair<fixed_tag, constant_t>,
      boost::mpl::pair<param_tag, parameter_t>
    > index_t_map;

public:
  typedef I index_tag;
  typedef typename boost::mpl::at<index_t_map, index_tag>::type index_t;
  typedef typename std::vector<index_t>::iterator iterator;
  typedef typename std::vector<index_t>::const_iterator const_iterator;

private:
  std::vector<index_t> indices;

public:
  ArrayIndex(const std::size_t _numIndices) : indices(_numIndices)
  {
    std::fill(indices.begin(), indices.end(), 0);
  }

  ArrayIndex(const std::size_t _numIndices, const constant_t* const _indices) :
    indices(_numIndices)
  {
    std::copy(_indices, _indices+_numIndices, indices.begin());
  }

  iterator begin()
  {
    return indices.begin();
  }

  const_iterator begin() const
  {
    return indices.begin();
  }

  iterator end()
  {
    return indices.end();
  }

  const_iterator end() const
  {
    return indices.end();
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
    return indices.size();
  }

  bool operator==(const ArrayIndex& i) const
  {
    return indices == i.indices;
  }

  bool operator<(const ArrayIndex& i) const
  {
    return indices < i.indices;
  }

  void swap(ArrayIndex& i)
  {
    std::swap(indices, i.indices);
  }
};

}

}

#endif
