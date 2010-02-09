#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_FUNCTION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_FUNCTION_HPP

#include <cstddef>
#include <set>
#include <boost/shared_ptr.hpp>
#include <simple_cfd/numeric/polynomial.hpp>
#include "array_fwd.hpp"
#include "scalar_reference.hpp"

namespace cfd
{

namespace detail
{

class TensorFunction
{
public:
  typedef boost::shared_ptr<TensorFunction> ref;
  typedef Polynomial<ScalarReference> polynomial_t;

  virtual std::size_t getTensorRank() const = 0;
  virtual std::size_t getTensorDimension() const = 0;
  virtual std::size_t numArrayIndices() const = 0;
  virtual std::size_t getArrayDimension(const std::size_t index) const = 0;
  virtual ref differentiate(const ScalarReference& reference) const = 0;
  virtual ~TensorFunction() {}
};

}

}
#endif
