#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_POLYNOMIAL_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_POLYNOMIAL_HPP

#include "tensor_array.hpp"
#include "tensor_array_table.hpp"
#include "scalar_placeholder.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace detail
{

class TensorArrayTablePolynomial : public TensorArrayTable<TensorArray::polynomial_t>
{
public:
  TensorArrayTablePolynomial(IndexGenerator& generator, const ArraySize& arraySize, 
    const TensorSize& tensorSize) : 
    TensorArrayTable<polynomial_t>(generator, arraySize, tensorSize)
  {
  }

  virtual TensorSize getTensorSize() const
  {
    return this->getTableTensorSize();
  }
  
  virtual TensorArrayRef derivative(const ScalarPlaceholder& x) const
  {
    if (x.isBound()) CFD_EXCEPTION("Cannot differentiate w.r.t. bound value.");

    TensorArrayRef result(*this);

  }
};

}

}

#endif
