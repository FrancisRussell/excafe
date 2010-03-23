#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_PLACEHOLDER_HPP

#include <vector>
#include <simple_cfd/exception.hpp>
#include "array_size.hpp"
#include "tensor_size.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayPlaceholder
{
private:
  long id;
  ArraySize arraySize;
  TensorSize tensorSize;

public:
  TensorArrayPlaceholder(const long _id, const ArraySize& _arraySize, const TensorSize& _tensorSize) :
    id(_id), arraySize(_arraySize), tensorSize(_tensorSize)
  {
  }

  bool operator<(const TensorArrayPlaceholder& p) const
  {
    if (id != p.id)
    {
      return id < p.id;
    }
    else
    {
      if (arraySize != p.arraySize || tensorSize != p.tensorSize)
        CFD_EXCEPTION("TensorArrayPlaceholder with same id but different sizes detected.");

      return false;
    }
  }

  bool operator==(const TensorArrayPlaceholder& p) const
  {
    if (id != p.id)
    {
      return false;
    }
    else
    {
      if (arraySize != p.arraySize || tensorSize != p.tensorSize)
        CFD_EXCEPTION("TensorArrayPlaceholder with same id but different sizes detected.");

      return true;
    }
  }

  long getID() const
  {
    return id;
  }

  ArraySize getArraySize() const
  {
    return arraySize;
  }

  TensorSize getTensorSize() const
  {
    return tensorSize;
  }
};

}

}
#endif
