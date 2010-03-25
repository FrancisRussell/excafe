#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_COLLECTIVE_SUM_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_COLLECTIVE_SUM_HPP

#include "tensor_size.hpp"
#include "tensor_array_collective.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayCollectiveSum : public TensorArrayCollective
{
public:
  TensorArrayCollectiveSum(const TensorSize& tensorSize) : 
    TensorArrayCollective(tensorSize)
  {
  }
};

}

}
#endif
