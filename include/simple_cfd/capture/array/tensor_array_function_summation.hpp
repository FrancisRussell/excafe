#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_SUMMATION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_SUMMATION_HPP

#include <cstddef>
#include <map>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>
#include "tensor_function.hpp"
#include "tensor_array_function_collective.hpp"
#include "tensor_array_function_helper.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "scalar_reference.hpp"
#include "index_incrementer.hpp"

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
      summand.setFunction(summand.getFunction()->differentiate(reference));
    }

    return TensorFunction::ref(new TensorArrayFunctionSummation(result)); 
  }

  virtual polynomial_t getPolynomial(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex) const
  {
    const std::map<ArrayIndexID, std::size_t> arrayIndexMap =
      TensorArrayFunctionHelper::indexToMap(this->arrayIndexParameters, arrayIndex);

    const std::map<TensorIndexID, std::size_t> tensorIndexMap =
      TensorArrayFunctionHelper::indexToMap(this->tensorIndexParameters, tensorIndex);

    const std::vector<std::size_t> unseenArrayExtentsVector(unseenArrayExtents.begin(), unseenArrayExtents.end());
    IndexIncrementer incrementer(unseenArrayIndexParameters, unseenTensorIndexParameters,
      unseenArrayExtentsVector, getTensorDimension());

    std::map<ArrayIndexID, std::size_t> unseenArrayIndex;
    std::map<TensorIndexID, std::size_t> unseenTensorIndex;

    incrementer.zero(unseenArrayIndex);
    incrementer.zero(unseenTensorIndex);

    polynomial_t poly(0.0);

    do
    {
      std::map<ArrayIndexID, std::size_t> fullArrayIndexMap(arrayIndexMap);
      fullArrayIndexMap.insert(unseenArrayIndex.begin(), unseenArrayIndex.end());

      std::map<TensorIndexID, std::size_t> fullTensorIndexMap(tensorIndexMap);
      fullTensorIndexMap.insert(unseenTensorIndex.begin(), unseenTensorIndex.end());

      BOOST_FOREACH(const call_t& call, this->operands)
      {
        poly += call.getFunction()->getPolynomial(
          TensorArrayFunctionHelper::getIndex(fullArrayIndexMap, call.getArrayIndex()),
          TensorArrayFunctionHelper::getIndex(fullTensorIndexMap, call.getTensorIndex()));
      }
    }
    while(!incrementer.increment(unseenArrayIndex, unseenTensorIndex));

    return poly;
  }

};

}

}

#endif
