#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_SUMMATION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_SUMMATION_HPP

#include <cstddef>
#include <map>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>
#include "tensor_function.hpp"
#include "tensor_array_function_collective.hpp"
#include "tensor_array_function_helper.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "scalar_reference.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionSummation : public TensorArrayFunctionCollective
{
public:
  TensorArrayFunctionSummation(const ArrayIndex<fixed_tag>& _arrayExtents, const std::size_t _rank, 
    const std::size_t _dimension) : TensorArrayFunctionCollective(_arrayExtents, _rank, _dimension)
  {
  }

  virtual TensorFunction::ref differentiate(const ScalarReference& reference)
  {
    TensorArrayFunctionSummation result(*this);
    BOOST_FOREACH(call_t& summand, result.operands)
    {
      boost::get<2>(summand) = boost::get<2>(summand)->differentiate(reference);
    }

    return TensorFunction::ref(new TensorArrayFunctionSummation(result)); 
  }

  virtual polynomial_t getPolynomial(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex) const
  {
    const std::map<ArrayIndexID, std::size_t> arrayIndexMap =
      TensorArrayFunctionHelper::indexToMap(this->arrayIndexParameters, arrayIndex);

    const std::map<TensorIndexID, std::size_t> tensorIndexMap =
      TensorArrayFunctionHelper::indexToMap(this->tensorIndexParameters, tensorIndex);

    polynomial_t poly(0.0);

    BOOST_FOREACH(const TensorArrayFunctionCollective::call_t& call, this->operands)
    {
      poly += boost::get<2>(call)->getPolynomial(
        TensorArrayFunctionHelper::getIndex(arrayIndexMap, boost::get<0>(call)),
        TensorArrayFunctionHelper::getIndex(tensorIndexMap, boost::get<1>(call)));
    }

    return poly;
  }

};

}

}

#endif
