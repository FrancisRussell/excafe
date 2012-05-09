#ifndef EXCAFE_NUMERIC_TENSOR_MATRIX_VIEW_HPP
#define EXCAFE_NUMERIC_TENSOR_MATRIX_VIEW_HPP

#include <cstddef>
#include <cassert>
#include "numeric_fwd.hpp"

namespace excafe
{

template<std::size_t D, typename T>
class TensorMatrixView
{
public:
  static const std::size_t  dimension = D;
  typedef T value_type;

private:
  Tensor<dimension, value_type>& tensor;

public:
  TensorMatrixView(Tensor<dimension, value_type>& _tensor) : tensor(_tensor)
  {
    assert(tensor.getRank() == 2);
  }

  value_type& operator()(const std::size_t row, const std::size_t col) 
  {
    TensorIndex index(TensorSize(2, dimension));
    index[0] = row;
    index[1] = col;

    return tensor[index];
  }

  const value_type& operator()(const std::size_t row, const std::size_t col) const
  {
    TensorIndex index(TensorSize(2, dimension));
    index[0] = row;
    index[1] = col;

    return tensor[index];
  }

  std::size_t getDimension() const
  {
    return tensor.getDimension();
  }

  std::size_t numRows() const
  {
    return getDimension();
  }

  std::size_t numCols() const
  {
    return getDimension();
  }
};

}

#endif
