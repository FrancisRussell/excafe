#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_HPP

#include "tensor_fwd.hpp"
#include "array_size.hpp"
#include "tensor_size.hpp"
#include <simple_cfd/numeric/polynomial.hpp>

namespace cfd
{

namespace detail
{

class TensorArray
{
public:
  typedef Polynomial<ScalarPlaceholder> polynomial_t;

  virtual TensorSize getTensorSize() const = 0;
  virtual TensorArrayRef derivative(const ScalarPlaceholder& x) const = 0;
  virtual ~TensorArray() {}
};

}

}

#endif
