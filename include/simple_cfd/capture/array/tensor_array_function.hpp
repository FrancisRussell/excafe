#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_HPP

#include <cstddef>
#include <set>
#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace detail
{

class TensorArrayFunction
{
public:
  typedef boost::shared_ptr<TensorArrayFunction> expr_ptr;

  virtual std::size_t getTensorRank() const = 0;
  virtual std::size_t getTensorDimension() const = 0;
  virtual std::size_t numArrayIndices() const = 0;
  virtual std::size_t getArrayDimension(const std::size_t index) const = 0;
  virtual std::size_t numScalarParameters() const = 0;
  virtual ~TensorArrayFunction() {}
};

}

}
#endif
