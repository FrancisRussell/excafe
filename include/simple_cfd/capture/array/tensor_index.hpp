#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_INDEX_HPP

#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>
#include <set>
#include <map>
#include <boost/mpl/map.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/at.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include "parameter_identifiers.hpp"
#include "array_traits.hpp"
#include "single_index.hpp"
#include <simple_cfd/exception.hpp>

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
  std::size_t dimension;
  std::vector<index_t> indices;

  TensorIndex(const std::size_t _rank, const std::size_t _dimension, const parameter_t* const _indices)
  {
    CFD_EXCEPTION("Cannot construct a constant array index from parameters.");
  }

public:
  TensorIndex(const std::size_t _rank, const std::size_t _dimension) : dimension(_dimension),
    indices(_rank)
  {
  }

  TensorIndex(const std::size_t _rank, const std::size_t _dimension, const constant_t* const _indices) :
    dimension(_dimension), indices(_rank)
  {
    std::copy(_indices, _indices+_rank, indices.begin());

    for(std::size_t i=0; i<_rank; ++i)
    {
      assert(_indices[i] < dimension);
    }
  }

  TensorIndex(const std::size_t _rank, const std::size_t _dimension, const TensorIndexID* const _indices)
  {
    CFD_EXCEPTION("Cannot construct a constant array index from tensor indices.");
  }

  TensorIndex(const TensorIndex<fixed_tag>& i) : dimension(i.getDimension()), indices(i.begin(), i.end())
  {
  }

  bool isParameterised() const
  {
    BOOST_STATIC_ASSERT(sizeof(index_tag) == 0);
    return false;
  }

  std::set<TensorIndexID> getReferencedParameters() const
  {
    BOOST_STATIC_ASSERT(sizeof(index_tag) == 0);
    return std::set<TensorIndexID>();
  }

  TensorIndex substituteLiterals(const std::map<TensorIndexID, std::size_t>& mapping) const
  {
    return *this;
  }

  TensorIndex head(const std::size_t n) const
  {
    assert(n <= indices.size());
    return TensorIndex(n, dimension, &indices[0]);
  }

  TensorIndex tail(const std::size_t n) const
  {
    assert(n <= indices.size());
    return TensorIndex(n, dimension, &indices[indices.size() - n]);
  }

  void prepend(const constant_t c)
  {
    indices.insert(indices.begin(), c);
  }

  void prepend(const TensorIndexID& id)
  {
    CFD_EXCEPTION("Cannot prepend a parameter to a constant array index.");
  }

  void append(const constant_t c)
  {
    indices.push_back(c);
  }

  void append(const TensorIndexID& id)
  {
    CFD_EXCEPTION("Cannot append a parameter to a constant array index.");
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
    return indices.size();
  }

  std::size_t getDimension() const
  {
    return dimension;
  }

  bool operator==(const TensorIndex& i) const
  {
    return dimension == i.dimension &&
           indices == i.indices;
  }

  bool operator<(const TensorIndex& i) const
  {
    return boost::make_tuple(dimension, indices) < boost::make_tuple(i.dimension, i.indices);
  }

  void swap(TensorIndex& i)
  {
    std::swap(dimension, i.dimension);
    std::swap(indices, i.indices);
  }
};

template<>
bool TensorIndex<param_tag>::isParameterised() const;

template<>
std::set<TensorIndexID> TensorIndex<param_tag>::getReferencedParameters() const;

template<>
TensorIndex<param_tag> TensorIndex<param_tag>::substituteLiterals(const std::map<TensorIndexID, std::size_t>& mapping) const;

template<>
TensorIndex<param_tag>::TensorIndex(const std::size_t _rank, const std::size_t dimension, const TensorIndexID* const _indices);

template<>
TensorIndex<param_tag>::TensorIndex(const std::size_t _rank, const std::size_t _dimension, const parameter_t* const _indices);

template<>
void TensorIndex<param_tag>::prepend(const TensorIndexID& id);

template<>
void TensorIndex<param_tag>::append(const TensorIndexID& id);

}

}

#endif
