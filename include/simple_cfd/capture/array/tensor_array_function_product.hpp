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

};


}

}

#endif
