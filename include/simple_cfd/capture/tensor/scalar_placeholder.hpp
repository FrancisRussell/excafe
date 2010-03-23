#ifndef SIMPLE_CFD_CAPTURE_TENSOR_SCALAR_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_SCALAR_PLACEHOLDER_HPP

#include "tensor_fwd.hpp"
#include "tensor_array_placeholder.hpp"
#include "index.hpp"
#include <simple_cfd/numeric/polynomial.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/operators.hpp>

namespace cfd
{

namespace detail
{

class ScalarPlaceholder
{
private:
  TensorArrayPlaceholder tensorArray;
  ArrayIndex arrayIndex;
  TensorIndex tensorIndex;

public:
  ScalarPlaceholder(const TensorArrayPlaceholder& _tensorArray, 
    const ArrayIndex& _arrayIndex, const TensorIndex& _tensorIndex) :
    tensorArray(_tensorArray), arrayIndex(_arrayIndex), tensorIndex(_tensorIndex)
  {
  }

  Polynomial<ScalarPlaceholder> operator-(const double d) const
  {
    return Polynomial<ScalarPlaceholder>(*this) - d;
  }

  Polynomial<ScalarPlaceholder> operator+(const double d) const
  {
    return Polynomial<ScalarPlaceholder>(*this) + d;
  }

  Polynomial<ScalarPlaceholder> operator*(const double d) const
  {
    return Polynomial<ScalarPlaceholder>(*this) * d;
  }

  Polynomial<ScalarPlaceholder> operator/(const double d) const
  {
    return Polynomial<ScalarPlaceholder>(*this) / d;
  }

  bool operator==(const ScalarPlaceholder& s) const
  {
    return tensorArray == s.tensorArray &&
           arrayIndex == s.arrayIndex &&
           tensorIndex == s.tensorIndex;
  }

  bool operator<(const ScalarPlaceholder& s) const
  {
    return boost::make_tuple(tensorArray, arrayIndex, tensorIndex) <
      boost::make_tuple(s.tensorArray, s.arrayIndex, s.tensorIndex);
  }
};

}

}

#endif
