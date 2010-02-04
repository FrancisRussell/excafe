#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_POLYNOMIAL_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_POLYNOMIAL_HPP

#include <cstddef>
#include <numeric>
#include <functional>
#include <vector>
#include "tensor_function.hpp"
#include "array_index.hpp"
#include "array_traits.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionPolynomial : public TensorFunction
{
private:
  typedef TensorFunction::polynomial_t polynomial_t;

  const ArrayIndex<fixed_tag> arrayExtents;
  const std::size_t rank;
  const std::size_t dimension;

  std::vector<polynomial_t> values;

  static std::size_t extent(const ArrayIndex<fixed_tag>& extents)
  {
    return std::accumulate(extents.begin(), extents.end(), 1, std::multiplies<ArrayIndex<fixed_tag>::constant_t>());
  }

  static std::size_t extent(const std::size_t rank, const std::size_t dimension)
  {
    std::size_t result=1;

    for(std::size_t i=0; i<rank; ++i)
      result *= dimension;

    return result;
  }

  std::size_t extent() const
  {
    return extent(arrayExtents) * extent(rank, dimension);
  }

public:
  TensorArrayFunctionPolynomial(const ArrayIndex<fixed_tag>& _arrayExtents, const std::size_t _rank, 
    const std::size_t _dimension) :
    arrayExtents(_arrayExtents), rank(_rank), dimension(_dimension), values(extent())
  {
  }
};

}

}

#endif
