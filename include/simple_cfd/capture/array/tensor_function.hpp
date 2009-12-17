#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_FUNCTION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_FUNCTION_HPP

#include <cstddef>
#include <cassert>
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/exception.hpp>
#include "array_expression.hpp"
#include "tensor_index.hpp"

namespace cfd
{

namespace detail
{

class TensorFunction : public ArrayExpression
{
private:
  const std::size_t rank;
  const std::size_t dimension;
  std::vector<Polynomial> functions;

  static std::size_t findTensorExtent(const std::size_t rank, const std::size_t dimension)
  {
    std::size_t extent = 1;

    for(std::size_t r=0; r<rank; ++r)
      extent *= dimension;

    return extent;
  }

public:
  TensorFunction(const std::size_t _rank, const std::size_t _dimension) :
    rank(_rank), dimension(_dimension), functions(findTensorExtent(rank, dimension))
  {
  }
  
  std::size_t getTensorRank() const
  {
    return rank;
  }

  std::size_t getTensorDimension() const
  {
    return dimension;
  }

  std::size_t getDimension(const std::size_t index) const
  {
    CFD_EXCEPTION("TensorFunctions have no array dimensions");
    return 0;
  }

  std::size_t numIndices() const
  {
    return 0;
  }

  Polynomial& operator[](const TensorIndex& index)
  {
    assert(rank == index.getRank());
    assert(dimension == index.getDimension());
    return functions[index.flatten()]; 
  }

  const Polynomial operator[](const TensorIndex& index) const
  {
    assert(rank == index.getRank());
    assert(dimension == index.getDimension());
    return functions[index.flatten()]; 
  }
};

}

}
#endif
