#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <simple_cfd/capture/tensor/tensor_placeholder.hpp>
#include <simple_cfd/capture/tensor/scalar_placeholder.hpp>

namespace cfd
{

namespace detail
{

TensorSize TensorPlaceholder::getTensorSize() const
{
  return tensorSize;
}

ScalarPlaceholder TensorPlaceholder::operator()(const TensorIndex::constant_t i) const
{
  const std::size_t rank = tensorSize.getRank();
  assert(rank == 1);

  TensorIndex tensorIndex(tensorSize);
  tensorIndex[0] = i;

  return ScalarPlaceholder(*this, tensorIndex);
}

ScalarPlaceholder TensorPlaceholder::operator[](const TensorIndex::constant_t i) const
{
  return (*this)(i);
}

bool TensorPlaceholder::operator==(const TensorPlaceholder& t) const
{
  return id == t.id && arrayIndices == t.arrayIndices && tensorSize == t.tensorSize;
}

bool TensorPlaceholder::operator<(const TensorPlaceholder& t) const
{
  return boost::make_tuple(id, arrayIndices, tensorSize) < 
    boost::make_tuple(t.id, t.arrayIndices, t.tensorSize);
}

}

}
