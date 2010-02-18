#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_SUMMATION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_SUMMATION_HPP

#include <cstddef>
#include <vector>
#include <map>
#include <utility>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include "tensor_function.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionSummation : public TensorFunction
{
private:
  typedef boost::tuple< ArrayIndex<param_tag>, TensorIndex<param_tag>, TensorFunction::ref > call_t;
  const ArrayIndex<fixed_tag> arrayExtents;
  const std::size_t rank;
  const std::size_t dimension;

  std::vector<ArrayIndexID> arrayIndexParameters;
  std::vector<TensorIndexID> tensorIndexParameters;
  std::vector<call_t> summands;

public:
  TensorArrayFunctionSummation(const ArrayIndex<fixed_tag>& _arrayExtents, const std::size_t _rank, 
    const std::size_t _dimension) : arrayExtents(_arrayExtents), rank(_rank), dimension(_dimension)
  {
    for(std::size_t i=0; i<arrayExtents.numIndices(); ++i)
       arrayIndexParameters.push_back(ArrayIndexID(i));

    for(std::size_t i=0; i<rank; ++i)
       tensorIndexParameters.push_back(TensorIndexID(i));
  }

  virtual std::size_t getTensorRank() const
  {
    return rank;
  }

  virtual std::size_t getTensorDimension() const
  {
    return dimension;
  }

  virtual std::size_t numArrayIndices() const
  {
    return arrayIndexParameters.size();
  }

  virtual std::size_t getArrayDimension(const std::size_t index) const
  {
    return arrayExtents[index];
  }

  virtual TensorFunction::ref differentiate(const ScalarReference& reference)
  {
    TensorArrayFunctionSummation result(*this);
    BOOST_FOREACH(call_t& summand, result.summands)
    {
      boost::get<2>(summand) = boost::get<2>(summand)->differentiate(reference);
    }

    return TensorFunction::ref(new TensorArrayFunctionSummation(result)); 
  }

  void addTerm(const ArrayIndex<param_tag>& arrayIndex, const TensorIndex<param_tag>& tensorIndex,
    const TensorFunction::ref function)
  {
    summands.push_back(boost::make_tuple(arrayIndex, tensorIndex, function));
  }

  ArrayIndex<param_tag> getDefaultArrayIndex() const
  {
    return ArrayIndex<param_tag>(arrayIndexParameters.size(), &arrayIndexParameters[0]);
  }

  TensorIndex<param_tag> getDefaultTensorIndex() const
  {
    return TensorIndex<param_tag>(tensorIndexParameters.size(), dimension, &tensorIndexParameters[0]);
  }
};

}

}

#endif
