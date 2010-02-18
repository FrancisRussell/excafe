#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_PRODUCT_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_PRODUCT_HPP

#include <cstddef>
#include <boost/tuple/tuple.hpp>
#include "tensor_function.hpp"
#include "tensor_array_function_collective.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionProduct : public TensorArrayFunctionCollective
{
public:
  TensorArrayFunctionProduct(const ArrayIndex<fixed_tag>& _arrayExtents, const std::size_t _rank, 
    const std::size_t _dimension) : TensorArrayFunctionCollective(_arrayExtents, _rank, _dimension)
  {
  }

  virtual TensorFunction::ref differentiate(const ScalarReference& reference)
  {
    //FIXME: implement me!
  }

  virtual polynomial_t getPolynomial(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex) const
  {
    const std::map<ArrayIndexID, std::size_t> arrayIndexMap =
      TensorArrayFunctionHelper::indexToMap(this->arrayIndexParameters, arrayIndex);

    const std::map<TensorIndexID, std::size_t> tensorIndexMap =
      TensorArrayFunctionHelper::indexToMap(this->tensorIndexParameters, tensorIndex);

    polynomial_t poly(1.0);

    BOOST_FOREACH(const TensorArrayFunctionCollective::call_t& call, this->operands)
    {
      poly *= boost::get<2>(call)->getPolynomial(
        TensorArrayFunctionHelper::getIndex(arrayIndexMap, boost::get<0>(call)),
        TensorArrayFunctionHelper::getIndex(tensorIndexMap, boost::get<1>(call)));
    }

    return poly;
  }
};


}

}

#endif
