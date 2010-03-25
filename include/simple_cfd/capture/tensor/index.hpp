#ifndef SIMPLE_CFD_CAPTURE_TENSOR_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_INDEX_HPP

#include <cstddef>
#include <deque>
#include <cassert>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/foreach.hpp>
#include "tensor_fwd.hpp"
#include "traits.hpp"
#include "array_size.hpp"
#include "tensor_size.hpp"
#include "index_expression.hpp"
#include "index_variable.hpp"

namespace cfd
{

namespace detail
{

template<typename T>
class Index
{
private:
  typedef T object_tag;
  typedef typename IndexProperties<object_tag>::size_type size_type;
  typedef typename IndexProperties<object_tag>::index_expression_t index_expression_t; 

  size_type size;
  std::deque<index_expression_t> indices;

  template<typename InputIterator>
  Index(const size_type& _size, const InputIterator begin, const InputIterator end) :
    size(_size), indices(begin, end)
  {
    assert(size.numIndices() == indices.size());
  }

public:
  typedef typename index_expression_t::constant_t constant_t;
  typedef typename index_expression_t::index_variable_t index_variable_t;

  static std::size_t flatten(const Index& index, const row_major_tag)
  {
    std::size_t flatIndex = 0;
    std::size_t multiplier = 1;

    for(int i=index.size.numIndices()-1; i>=0; --i)
    {
      flatIndex += multiplier * index.indices[i].toConstant();
      multiplier *= index.size.getLimit(i);
    }

    return flatIndex;
  }
  
  static Index unflatten(const size_type& size, const std::size_t flatIndex, const row_major_tag)
  {
    Index index(size);
    std::size_t multiplier = size.getExtent();
    std::size_t remainingIndex = flatIndex;

    for(std::size_t i=0; i<size.numIndices(); ++i)
    {
      multiplier /= size.getLimit(i);
      index[i] = remainingIndex / multiplier;
      remainingIndex %= multiplier;
    }

    return index;
  }

  Index(const size_type& _size) : size(_size), indices(size.numIndices())
  {
  }

  index_expression_t& operator[](const std::size_t index)
  {
    assert(index < indices.size());
    return indices[index];
  }

  const index_expression_t& operator[](const std::size_t index) const
  {
    assert(index < indices.size());
    return indices[index];
  }

  size_type getSize() const
  {
    return size;
  }

  bool operator==(const Index& i) const
  {
    return size == i.size &&
           indices == i.indices;
  }

  bool operator<(const Index& i) const
  {
    return boost::make_tuple(size, indices) < boost::make_tuple(i.size, i.indices);
  }

  bool allConstant() const
  {
    bool allConstant = true;
    BOOST_FOREACH(const index_expression_t& e, indices)
    {
      allConstant = allConstant && e.isConstant();
    }
    return allConstant;
  }

  bool allVariable() const
  {
    bool allVariable = true;
    BOOST_FOREACH(const index_expression_t& e, indices)
    {
      allVariable = allVariable && e.isVariable();
    }
    return allVariable;
  }

  Index prepend(const index_variable_t& v) const
  {
    Index result(*this);
    result.size = result.size.prepend(v.getLimit());
    result.indices.push_front(v);
    return result;
  }

  Index append(const index_variable_t& v)
  {
    Index result(*this);
    result.size = result.size.append(v.getLimit());
    result.indices.push_back(v);
    return result;
  }

  Index head(const std::size_t n) const
  {
    assert(n <= size.numIndices());
    return Index(size.head(n), indices.begin(), indices.begin()+n);
  }

  Index tail(const std::size_t n) const
  {
    assert(n <= size.numIndices());
    return Index(size.tail(n), indices.end()-n, indices.end());
  }
};

}

}

#endif
