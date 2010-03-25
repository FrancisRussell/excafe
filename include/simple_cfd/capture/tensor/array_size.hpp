#ifndef SIMPLE_CFD_CAPTURE_TENSOR_ARRAY_SIZE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_ARRAY_SIZE_HPP

#include <cstddef>
#include <cassert>
#include <deque>
#include <numeric>
#include <boost/operators.hpp>
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace detail
{

class ArraySize : boost::equality_comparable<ArraySize>
{
private:
  typedef std::size_t value_type;
  typedef std::deque<value_type>::iterator iterator;
  typedef std::deque<value_type>::const_iterator const_iterator;

  std::deque<value_type> limits;

  template<typename InputIterator>
  ArraySize(const InputIterator begin, const InputIterator end) : limits(begin, end)
  {
  }

public:
  ArraySize() : limits(0)
  {
  }

  ArraySize(const std::size_t _indexCount) : limits(_indexCount)
  {
  }

  ArraySize(const std::size_t _indexCount, const value_type* const _limits) : limits(_indexCount)
  {
    std::copy(_limits, _limits+_indexCount, limits.begin());
  }

  std::size_t numIndices() const
  {
    return limits.size();
  }
  
  std::size_t getExtent() const
  {
    const std::size_t extent = std::accumulate(begin(), end(), 1, std::multiplies<std::size_t>());
    return extent;
  }

  value_type& operator[](const std::size_t i)
  {
    assert(i<limits.size());
    return limits[i];
  }

  const value_type& operator[](const std::size_t i) const
  {
    assert(i<limits.size());
    return limits[i];
  }

  std::size_t getLimit(const std::size_t index) const
  {
    return (*this)[index];
  }

  iterator begin()
  {
    return limits.begin();
  }

  const_iterator begin() const
  {
    return limits.begin();
  }

  iterator end()
  {
    return limits.end();
  }

  const_iterator end() const
  {
    return limits.end();
  }

  bool operator==(const ArraySize& s) const
  {
    return limits == s.limits;
  }

  bool operator<(const ArraySize& s) const
  {
    return limits < s.limits;
  }

  ArraySize prepend(const std::size_t d) const
  {
    ArraySize result(*this);
    result.limits.push_front(d);
    return result;
  }

  ArraySize append(const std::size_t d) const
  {
    ArraySize result(*this);
    result.limits.push_back(d);
    return result;
  }

  ArraySize head(const std::size_t n) const
  {
    assert(n <= limits.size());
    return ArraySize(limits.begin(), limits.begin()+n);
  }

  ArraySize tail(const std::size_t n) const
  {
    assert(n <= limits.size());
    return ArraySize(limits.end()-n, limits.end());
  }

};

}

}

#endif
