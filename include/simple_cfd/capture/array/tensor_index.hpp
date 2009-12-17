#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_INDEX_HPP

#include <cstddef>
#include <vector>
#include <cassert>
#include <algorithm>

namespace cfd
{

namespace detail
{

class TensorIndex
{
private:
  std::size_t rank;
  std::size_t dimension;
  std::vector<std::size_t> indices;

public:
  TensorIndex(const std::size_t _rank, const std::size_t _dimension) : rank(_rank), dimension(_dimension),
    indices(rank)
  {
  }

  TensorIndex(const std::size_t _rank, const std::size_t _dimension, const std::size_t* const _indices) :
    rank(_rank), dimension(_dimension), indices(rank)
  {
    std::copy(_indices, _indices+rank, indices.begin());

    for(std::size_t i=0; i<indices.size(); ++i)
    {
      assert(indices[i] < dimension);
    }
  }

  std::size_t& operator[](const std::size_t index)
  {
    assert(index < indices.size());
    return indices[index];
  }

  std::size_t operator[](const std::size_t index) const
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

  std::size_t flatten() const
  {
    std::size_t index = 0;

    for(std::size_t i=0; i<indices.size(); ++i)
    {
      index *= dimension;
      index += indices[i];
    }

    return index;
  }

  bool operator==(const TensorIndex& i) const
  {
    return rank == i.rank &&
           dimension == i.dimension &&
           indices == i.indices;
  }

  bool operator<(const TensorIndex& i) const
  {
    if (rank < i.rank) return true;
    if (rank == i.rank && dimension < i.dimension) return true;
    if (rank == i.rank && dimension == i.dimension && indices < i.indices) return true;
    return false;
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
