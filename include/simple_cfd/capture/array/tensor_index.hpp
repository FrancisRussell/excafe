#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_INDEX_HPP

#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>
#include <boost/mpl/map.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/at.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/foreach.hpp>
#include "parameter_identifiers.hpp"
#include "array_traits.hpp"

namespace cfd
{

namespace detail
{

template<typename I>
class TensorIndex
{
public:
  typedef std::size_t constant_t;
  typedef SingleIndex<TensorIndexID> parameter_t;

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
  std::size_t rank;
  std::size_t dimension;
  std::vector<index_t> indices;

public:
  TensorIndex(const std::size_t _rank, const std::size_t _dimension) : rank(_rank), dimension(_dimension),
    indices(rank)
  {
    std::fill(indices.begin(), indices.end(), 0);
  }

  TensorIndex(const std::size_t _rank, const std::size_t _dimension, const constant_t* const _indices) :
    rank(_rank), dimension(_dimension), indices(rank)
  {
    std::copy(_indices, _indices+rank, indices.begin());

    for(std::size_t i=0; i<rank; ++i)
    {
      assert(_indices[i] < dimension);
    }
  }

  TensorIndex(const std::size_t _rank, const std::size_t _dimension, const parameter_t* const _indices) :
    rank(_rank), dimension(_dimension), indices(rank)
  {
    std::copy(_indices, _indices+rank, indices.begin());
  }

  bool isParameterised() const
  {
    bool parameterised = false;

    BOOST_FOREACH(const index_t& index, indices)
    {
      parameterised |= index.isParameter();
    }
    
    return parameterised;
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

  std::size_t getRank() const
  {
    return rank;
  }

  std::size_t getDimension() const
  {
    return dimension;
  }

  bool operator==(const TensorIndex& i) const
  {
    return rank == i.rank &&
           dimension == i.dimension &&
           indices == i.indices;
  }

  bool operator<(const TensorIndex& i) const
  {
    return boost::make_tuple(rank, dimension, indices) < boost::make_tuple(i.rank, i.dimension, i.indices);
  }

  void swap(TensorIndex& i)
  {
    std::swap(rank, i.rank);
    std::swap(dimension, i.dimension);
    std::swap(indices, i.indices);
  }
};

}

}

#endif
