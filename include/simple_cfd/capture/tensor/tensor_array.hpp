#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_HPP

#include "array_size.hpp"
#include "tensor_size.hpp"
#include "scalar_placeholder.hpp"
#include <simple_cfd/numeric/polynomial.hpp>

namespace cfd
{

namespace detail
{

class TensorArray
{
public:
  typedef Polynomial<ScalarPlaceholder> polynomial_t;

  virtual ArraySize getArraySize() const = 0;
  virtual TensorSize getTensorSize() const = 0;
  virtual ~TensorArray() {}
};

}

}

#endif
