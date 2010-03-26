#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_COLLECTIVE_PRODUCT_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_COLLECTIVE_PRODUCT_HPP

#include <cstddef>
#include "tensor_fwd.hpp"
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

  TensorArrayCollectiveProduct(const TensorArrayCollectiveProduct& t) :
    TensorArrayCollective(t)
  {
  }

  virtual TensorArrayRef derivative(const ScalarPlaceholder& x) const
  {
    TensorArrayCollectiveSum result(this->getTensorSize()); 

    for(std::size_t i=0; i<numOperands(); ++i)
    {
      TensorArrayCollectiveProduct subTerm(*this);
      IndexedTensor& differentiated = subTerm.operand(i);
      differentiated = IndexedTensor(differentiated.getTensorRef()->derivative(x),
        differentiated.getIndex());
      result.addOperand(result.getVisibleIndices(), subTerm);
    }

    return result;
  }

};

}

}
#endif
