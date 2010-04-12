#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_POLYNOMIAL_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_TABLE_POLYNOMIAL_HPP

#include <boost/foreach.hpp>
#include "tensor_array.hpp"
#include "tensor_array_table.hpp"
#include "scalar_placeholder.hpp"

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

  TensorArrayTablePolynomial(IndexGenerator& generator, const ArrayIndex& _arrayIndices, 
    const TensorSize& _tensorSize) :
    TensorArrayTable<polynomial_t>(generator, _arrayIndices, _tensorSize)
  {
  }

  virtual TensorArrayRef derivative(const ScalarPlaceholder& x) const
  {
    TensorArrayTablePolynomial result(*this);
    BOOST_FOREACH(TensorArray::polynomial_t& element, result)
    {
      element = element.derivative(x);
    }

    return TensorArrayRef::cloneFrom(result);
  }
};

}

}

#endif
