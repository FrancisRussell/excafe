#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_VISITOR_HPP

#include <map>
#include <cstddef>
#include "array_fwd.hpp"
#include "parameter_identifiers.hpp"
#include "tensor_function.hpp"

namespace cfd
{

namespace detail
{

template<typename E>
class TensorArrayFunctionVisitor
{
public:
  typedef E element_t;

  virtual void visit(const TensorArrayFunction<element_t>& parent,
    const std::map<ArrayIndexID, std::size_t>& arrayIndex, 
    const std::map<TensorIndexID, std::size_t>& tensorIndex,
    element_t& value) = 0;
  
  virtual ~TensorArrayFunctionVisitor() {}
};

}

}

#endif
