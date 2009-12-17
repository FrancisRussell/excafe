#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_FUNCTION_ARRAY_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_FUNCTION_ARRAY_HPP

namespace cfd
{

namespace detail
{

class TensorFunctionArray : public ArrayExpression
{
public:
  ArrayIndex arrayDimensions;
  const std::size_t rank;
  const std::size_t dimension;

public:
  TensorFunctionArray(const ArrayIndex& i, const std::size_t _rank, const std::size_t _dimension) :
    arrayDimensions(i), rank(_rank), dimension(_dimension)
  {
  }

  virtual std::size_t getTensorRank() const
  {
    return rank;
  }

  virtual std::size_t getTensorDimension() const
  {
    return dimension;
  }

  virtual std::size_t numArrayIndices() const
  {
    return arrayDimensions.getNumIndices();
  }

  virtual std::size_t getArrayDimension(const std::size_t index) const
  {
    return arrayDimensions.getArrayDimension(index);
  }
};

}

}

#endif
