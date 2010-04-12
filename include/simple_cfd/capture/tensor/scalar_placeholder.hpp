#ifndef SIMPLE_CFD_CAPTURE_TENSOR_SCALAR_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_SCALAR_PLACEHOLDER_HPP

#include <cassert>
#include <utility>
#include <boost/operators.hpp>
#include "tensor_fwd.hpp"
#include "tensor_array_ref.hpp"
#include "index.hpp"
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace detail
{

class ScalarPlaceholder : boost::equality_comparable<ScalarPlaceholder>
{
private:
  TensorArrayRef tensor;
  TensorIndex tensorIndex;

public:
  ScalarPlaceholder(const TensorArrayRef& _tensor, const TensorIndex& _tensorIndex) :
    tensor(_tensor), tensorIndex(_tensorIndex)
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
    return tensor == s.tensor &&
           tensorIndex == s.tensorIndex;
  }

  bool operator<(const ScalarPlaceholder& s) const
  {
    return std::make_pair(tensor, tensorIndex) <
      std::make_pair(s.tensor, s.tensorIndex);
  }

  TensorArrayRef getTensorArrayRef() const
  {
    return tensor;
  }

  TensorIndex getTensorIndex() const
  {
    return tensorIndex;
  }
};

}

}

#endif
