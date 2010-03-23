#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_SIZE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_SIZE_HPP

#include <cstddef>
#include <cassert>
#include <utility>
#include <boost/operators.hpp>

namespace cfd
{

namespace detail
{

class TensorSize : boost::equality_comparable<TensorSize>
{
private:
  std::size_t rank;
  std::size_t dimension;

public:
  TensorSize(const std::size_t _rank, const std::size_t _dimension) :
    rank(_rank), dimension(_dimension)
  {
  }

  std::size_t getRank() const
  {
    return rank;
  }

  std::size_t getDimension() const
  {
    return dimension;
  }

  std::size_t numIndices() const
  {
    return getRank();
  }

  std::size_t getLimit(const std::size_t index) const
  {
    assert(index < getRank());
    return getDimension();
  }

  std::size_t getExtent() const
  {
    std::size_t extent = 1;

    for(std::size_t r=0; r<rank; ++r)
      extent *= dimension;

    return extent;
  }

  bool operator==(const TensorSize& s) const
  {
    return rank == s.rank && dimension == s.dimension;
  }

  bool operator<(const TensorSize& s) const
  {
    return std::make_pair(rank,dimension) < std::make_pair(s.rank, s.dimension);
  }
};

}

}

#endif
