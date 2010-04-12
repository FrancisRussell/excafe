#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_MATRIX_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_MATRIX_HPP

#include "tensor_array.hpp"
#include "tensor_array_ref.hpp"
#include "index_generator.hpp"
#include <cstddef>
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace detail
{

class TensorArrayMatrixRight
{
private:
  TensorArrayRef ref;

  std::size_t getDimension() const
  {
    return ref->getTensorSize().getDimension();
  }

public:
  TensorArrayMatrixRight(const TensorArrayRef& _ref) : ref(_ref)
  {
    if (ref->getTensorSize().getRank() != 2) CFD_EXCEPTION("Can only wrap tensor of rank 2 as a matrix.");
  }

  std::size_t numRows() const
  {
    return getDimension();
  }

  std::size_t numCols() const
  {
    return getDimension();
  }

  const ScalarPlaceholder operator()(const std::size_t row, const std::size_t col) const
  {
    TensorIndex index(ref->getTensorSize());
    index[0] = row;
    index[1] = col;
    return ScalarPlaceholder(ref, index);
  }
};

class TensorArrayMatrixLeft
{
private:
  IndexGenerator generator;
  ArraySize arraySize;
  TensorSize tensorSize;
  TensorArrayTablePolynomial matrix;

public:
  TensorArrayMatrixLeft(const std::size_t n) : arraySize(0), tensorSize(2, n), 
    matrix(generator, arraySize, tensorSize)
  {
  }

  std::size_t numRows() const
  {
    return tensorSize.getDimension();
  }

  std::size_t numCols() const
  {
    return tensorSize.getDimension();
  }

  TensorArray::polynomial_t& operator()(const std::size_t row, const std::size_t col)
  {
    const ArrayIndex arrayIndex(arraySize);
    TensorIndex tensorIndex(tensorSize);
    tensorIndex[0] = row;
    tensorIndex[1] = col;

    return matrix(arrayIndex, tensorIndex);
  }

  TensorArrayRef toTensorArrayRef() const
  {
    return TensorArrayRef::cloneFrom(matrix);
  }
};

}

}

#endif
