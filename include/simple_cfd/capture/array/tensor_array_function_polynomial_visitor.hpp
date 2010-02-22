#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_POLYNOMIAL_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_POLYNOMIAL_VISITOR_HPP

#include <map>
#include <cstddef>
#include "parameter_identifiers.hpp"
#include "tensor_function.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionPolynomialVisitor
{
public:
  virtual void visit(const TensorArrayFunctionPolynomial& parent,
    const std::map<ArrayIndexID, std::size_t>& arrayIndex, 
    const std::map<TensorIndexID, std::size_t>& tensorIndex,
    TensorFunction::polynomial_t& value) = 0;
  
  virtual ~TensorArrayFunctionPolynomialVisitor() {}
};

}

}

#endif
