#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_OPERATORS_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_OPERATORS_HPP

namespace cfd
{

namespace detail
{

template<typename S, typename T, typename U>
struct ImmutableArithmetic
{
  friend S operator+(const T& t, const U& u) { return S(t) + u; }
  friend S operator-(const T& t, const U& u) { return S(t) - u; }
  friend S operator/(const T& t, const U& u) { return S(t) / u; }
  friend S operator*(const T& t, const U& u) { return S(t) * u; }
};

}

}

#endif
