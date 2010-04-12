#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_PLACEHOLDER_HPP

#include <cstddef>
#include "index.hpp"
#include "index_generator.hpp"
#include "array_size.hpp"
#include "tensor_size.hpp"
#include "tensor_array.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace detail
{

class TensorArrayPlaceholder : public TensorArray
{
private:
  TensorSize tensorSize;

public:
  TensorArrayPlaceholder(const TensorSize& _tensorSize) :
    tensorSize(_tensorSize)
  {
  }

  TensorSize getTensorSize() const
  {
    return tensorSize;
  }

  TensorArrayRef derivative(const ScalarPlaceholder& x) const
  {
    const ArraySize nullArraySize(0);
    const ArrayIndex nullArrayIndex(nullArraySize);
    IndexGenerator g;
    TensorArrayTablePolynomial result(g, nullArraySize, tensorSize);

    if (x.getTensorArrayRef().get() == this)
    {
      if (!x.getTensorIndex().allConstant()) 
        CFD_EXCEPTION("Can only differentiate w.r.t a tensor indexed with a constant value.");

      result(nullArrayIndex, x.getTensorIndex()) = TensorArray::polynomial_t(1.0);
    }

    return TensorArrayRef::cloneFrom(result);
  }

};

}

}
#endif
