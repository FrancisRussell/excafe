#ifndef SIMPLE_CFD_CAPTURE_TENSOR_INDEXED_TENSOR_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_INDEXED_TENSOR_HPP

#include <utility>
#include <cassert>
#include "tensor_fwd.hpp"
#include "tensor_array_ref.hpp"
#include "index.hpp"

namespace cfd
{

namespace detail
{

class IndexedTensor
{
private:
  TensorArrayRef tensor;
  TensorIndex tensorIndex;

public:
  IndexedTensor(const TensorArrayRef& _tensor, const TensorIndex& _tensorIndex) :
    tensor(_tensor), tensorIndex(_tensorIndex)
  {
    assert(tensor->getTensorSize() == tensorIndex.getSize());
  }

  bool operator==(const IndexedTensor& s) const
  {
    return tensor == s.tensor &&
           tensorIndex == s.tensorIndex;
  }

  bool operator<(const IndexedTensor& s) const
  {
    return std::make_pair(tensor, tensorIndex) <
      std::make_pair(s.tensor, s.tensorIndex);
  }
};

}

}

#endif
