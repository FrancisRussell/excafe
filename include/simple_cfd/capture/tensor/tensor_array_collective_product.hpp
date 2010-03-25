#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_COLLECTIVE_PRODUCT_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_COLLECTIVE_PRODUCT_HPP

#include "tensor_size.hpp"
#include "tensor_array_collective.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayCollectiveProduct : public TensorArrayCollective
{
public:
  TensorArrayCollectiveProduct(const TensorSize& tensorSize) : 
    TensorArrayCollective(tensorSize)
  {
  }
};

}

}
#endif
