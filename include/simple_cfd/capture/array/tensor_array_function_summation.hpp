#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_SUMMATION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_SUMMATION_HPP

#include <cstddef>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>
#include "tensor_function.hpp"
#include "tensor_array_function_collective.hpp"
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
};

}

}

#endif
