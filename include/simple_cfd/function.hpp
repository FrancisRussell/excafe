#ifndef SIMPLE_CFD_FUNCTION_HPP
#define SIMPLE_CFD_FUNCTION_HPP

#include "simple_cfd_fwd.hpp"

namespace cfd
{

template<unsigned int D, unsigned int R>
class Function
{
public:
  static const int dimension = D;
  static const int rank = R;
  typedef typename Tensor<D, R>::value_type value_type;

  virtual Tensor<D, R> evaluate(const vertex<dimension>& v) const = 0;
};

}

#endif
