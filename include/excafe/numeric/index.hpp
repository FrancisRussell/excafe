#ifndef EXCAFE_CAPTURE_ASSEMBLY_INDEX_HPP
#define EXCAFE_CAPTURE_ASSEMBLY_INDEX_HPP

#include <cstddef>
#include <deque>
#include <cassert>
#include <ostream>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/foreach.hpp>
#include "numeric_fwd.hpp"
#include "traits.hpp"
#include "tensor_size.hpp"

namespace excafe
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
  static std::size_t flatten(const Index& index, const row_major_tag)
  {
    std::size_t flatIndex = 0;
    std::size_t multiplier = 1;

    for(int i=index.size.numIndices()-1; i>=0; --i)
    {
      flatIndex += multiplier * index.indices[i];
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

  std::size_t numIndices() const
  {
    assert(size.numIndices() == indices.size());
    return size.numIndices();
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

  void write(std::ostream& o) const
  {
    o << "(";

    for(std::size_t i=0; i<size.numIndices(); ++i)
    {
      o << indices[i];

      if (i != size.numIndices()-1)
        o << ", ";
    }

    o << ")";
  }
};

template<typename T>
std::ostream& operator<<(std::ostream& o, const Index<T>& index) 
{
  index.write(o);
  return o;
}

}

}

#endif
