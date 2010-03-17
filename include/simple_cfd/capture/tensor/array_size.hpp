#ifndef SIMPLE_CFD_CAPTURE_TENSOR_ARRAY_SIZE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_ARRAY_SIZE_HPP

#include <cstddef>
#include <cassert>
#include <vector>

namespace cfd
{

namespace detail
{

class ArraySize
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

  ArraySize(const std::size_t _indexCount) : limits(indexCount)
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
};

}

}

#endif
