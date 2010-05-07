#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_OPERATORS_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_OPERATORS_HPP

#include "assembly_fwd.hpp"
#include <simple_cfd/numeric/polynomial.hpp>

namespace cfd
{

namespace detail
{

template<typename T>
class ScalarPlaceholderOperators
{
private:
  typedef T child_t;
  typedef Polynomial<ScalarPlaceholder> polynomial_t;

  const child_t& toChild() const
  {
    return *static_cast<const child_t*>(this);
  }

  const polynomial_t toPolynomial() const
  {
    return polynomial_t(ScalarPlaceholder(toChild()));
  }

public:
  polynomial_t operator-(const double d) const
  {
    return toPolynomial() - d;
  }
};

}

}

#endif
