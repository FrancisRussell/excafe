#ifndef SIMPLE_CFD_CAPTURE_TENSOR_ARRAY_SIZE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_ARRAY_SIZE_HPP

#include <cstddef>
#include <cassert>
#include <vector>
#include <numeric>
#include <boost/operators.hpp>

namespace cfd
{

namespace detail
{

class ArraySize : boost::equality_comparable<ArraySize>
{
private:
  typedef std::size_t value_type;
  typedef std::vector<value_type>::iterator iterator;
  typedef std::vector<value_type>::const_iterator const_iterator;

  std::vector<value_type> limits;

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
    std::accumulate(begin(), end(), 1, std::multiplies<std::size_t>());
    std::size_t extent = 1;

    for(std::size_t i=0; i<limits.size(); ++i)
      extent *= limits[i];

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
};

}

}

#endif
