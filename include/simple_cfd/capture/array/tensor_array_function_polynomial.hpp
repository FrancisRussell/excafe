#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_POLYNOMIAL_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_POLYNOMIAL_HPP

#include "tensor_function.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionPolynomial : public TensorFunction
{
private:
  const ArrayIndex arrayExtents;
  const std::size_t rank;
  const std::size_t dimension;

public:
  TensorArrayFunctionPolynomial(const ArrayIndex& _arrayExtents, const std::size_t _rank, const std::size_t _dimension)
  {
  }
};

}

}

#endif
