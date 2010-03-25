#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_SIZE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_SIZE_HPP

#include <cstddef>
#include <cassert>
#include <utility>
#include <boost/operators.hpp>
#include <simple_cfd/exception.hpp>


namespace cfd
{

namespace detail
{

class TensorSize : boost::equality_comparable<TensorSize>
{
private:
  std::size_t rank;
  std::size_t dimension;

  TensorSize push(const std::size_t d) const
  {
    if (d != dimension) CFD_EXCEPTION("Can only increase tensor sizes by values of dimension."); 
    TensorSize result(*this);
    ++result.rank;
    return result;
  }

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

  TensorSize prepend(const std::size_t d) const
  {
    return push(d);
  }

  TensorSize append(const std::size_t d) const
  {
    return push(d);
  }

  TensorSize head(const std::size_t n) const
  {
    assert(n <= rank);
    return TensorSize(rank-n, dimension);
  }

  TensorSize tail(const std::size_t n) const
  {
    assert(n <= rank);
    return TensorSize(rank-n, dimension);
  }
};

}

}

#endif
