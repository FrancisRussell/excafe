#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_COLLECTIVE_SUM_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_COLLECTIVE_SUM_HPP

#include <boost/foreach.hpp>
#include "tensor_fwd.hpp"
#include "tensor_size.hpp"
#include "indexed_tensor.hpp"
#include "tensor_array_collective.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayCollectiveSum : public TensorArrayCollective
{
public:
  TensorArrayCollectiveSum(const TensorArrayCollectiveSum& t) :
    TensorArrayCollective(t)
  {
  }

  TensorArrayCollectiveSum(const TensorSize& tensorSize) : 
    TensorArrayCollective(tensorSize)
  {
  }

  virtual TensorArrayRef derivative(const ScalarPlaceholder& x) const
  {
    TensorArrayCollectiveSum result(*this);

    BOOST_FOREACH(IndexedTensor& indexed, result)
    {
      indexed = IndexedTensor(indexed.getTensorRef()->derivative(x), indexed.getIndex());
    }
    return TensorArrayRef::cloneFrom(result);
  }
};

}

}
#endif
