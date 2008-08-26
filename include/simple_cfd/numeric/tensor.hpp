#ifndef SIMPLE_CFD_NUMERIC_TENSOR_HPP
#define SIMPLE_CFD_NUMERIC_TENSOR_HPP

namespace cfd
{

namespace detail
{

template<unsigned  X, unsigned Y>
struct Power
{
  static const unsigned int value = Power<X*Y, Y-1>::value;
};

template<unsigned X>
class Power<X, 0>
{
  static const unsigned int value = 1;
};

}

template<unsigned int D, unsigned int R, typename T>
class Tensor
{
public:
  static const unsigned int dimension = D;
  static const unsigned int rank = R;
  static const unsigned int size = detail::Power<D, R>::value;
};

}

#endif
