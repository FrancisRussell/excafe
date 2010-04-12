#include <cassert>
#include <cstddef>
#include <simple_cfd/capture/tensor/index.hpp>
#include <simple_cfd/capture/tensor/scalar_placeholder.hpp>
#include <simple_cfd/capture/tensor/tensor_array_ref.hpp>

namespace cfd
{

namespace detail
{

ScalarPlaceholder TensorArrayRef::operator[](const TensorIndex::constant_t i) const
{
  const TensorSize tensorSize((*this)->getTensorSize()); 
  const std::size_t rank = tensorSize.getRank();
  assert(rank == 1);

  TensorIndex tensorIndex(tensorSize);
  tensorIndex[0] = i;

  return ScalarPlaceholder(*this, tensorIndex);
}

}

}
