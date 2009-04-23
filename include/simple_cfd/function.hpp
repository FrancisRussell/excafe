#ifndef SIMPLE_CFD_FUNCTION_HPP
#define SIMPLE_CFD_FUNCTION_HPP

#include "simple_cfd_fwd.hpp"

namespace cfd
{

template<unsigned int D, unsigned int R, typename T>
class Function
{
public:
  static const int dimension = D;
  static const int rank = R;
  typedef T value_type;

  virtual Tensor<D, R, T> evaluate(const vertex<dimension>& v) const = 0;
};

}

#endif
