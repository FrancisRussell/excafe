#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_REFERENCES_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_REFERENCES_HPP

#include <cstddef>
#include <functional>
#include <algorithm>
#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <iterator>
#include <boost/foreach.hpp>
#include <simple_cfd/exception.hpp>
#include "tensor_function.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "array_traits.hpp"
#include "tensor_array_function.hpp"
#include "tensor_array_function_visitor.hpp"
#include "tensor_array_function_helper.hpp"
#include "tensor_array_function_call.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionReferences : public TensorArrayFunction<TensorArrayFunctionCall>
{
private:
  typedef TensorArrayFunctionCall element_t;
  typedef TensorArrayFunction<element_t> parent_t;
  typedef TensorArrayFunctionVisitor<element_t> visitor_t;

  class Differentiator : public visitor_t
  {
  private:
    const ScalarReference variable;

  public:
    Differentiator(const ScalarReference& _variable) : variable(_variable)
    {
      assert(!variable.isBound() && !variable.isParameterised());
    }

    void visit(const TensorArrayFunction<element_t>& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndexMap,
      const std::map<TensorIndexID, std::size_t>& tensorIndexMap,
      element_t& value)
    {
      value.setFunction(value.getFunction()->differentiate(variable));
    }
  };

public:
  TensorArrayFunctionReferences(const ArrayIndex<fixed_tag>& _arrayExtents, const std::size_t _rank, 
    const std::size_t _dimension) : parent_t(_arrayExtents, _rank, _dimension)
  {
  }

  virtual polynomial_t getPolynomial(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex) const
  {
    const std::map<ArrayIndexID, std::size_t> arrayIndexMap =
      TensorArrayFunctionHelper::indexToMap(arrayIndexParameters, arrayIndex);

    const std::map<TensorIndexID, std::size_t> tensorIndexMap =
      TensorArrayFunctionHelper::indexToMap(tensorIndexParameters, tensorIndex);

    const element_t call = (*this)(arrayIndexMap, tensorIndexMap);

    //TODO: assert that this polynomial isn't parameterised

    return call.getFunction()->getPolynomial(
      TensorArrayFunctionHelper::getIndex(arrayIndexMap, call.getArrayIndex()),
      TensorArrayFunctionHelper::getIndex(tensorIndexMap, call.getTensorIndex()));
  }

  virtual TensorFunction::ref differentiate(const ScalarReference& reference)
  {
    if (reference.isBound()) 
        CFD_EXCEPTION("Cannot differentiate with respect to bound reference.");

    if (reference.isParameterised()) 
      CFD_EXCEPTION("Cannot differentiate with respect to parameterised reference.");
    
    TensorArrayFunctionReferences result(*this);

    Differentiator differentiator(reference);
    result.accept(differentiator);

    return TensorFunction::ref(new TensorArrayFunctionReferences(result));
  }
};

}

}

#endif
