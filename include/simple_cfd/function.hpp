#ifndef SIMPLE_CFD_FUNCTION_HPP
#define SIMPLE_CFD_FUNCTION_HPP

#include "simple_cfd_fwd.hpp"

namespace cfd
{

template<unsigned int D>
class Function
{
public:
  static const int dimension = D;
  typedef typename Tensor<dimension>::value_type value_type;

  virtual Tensor<dimension> evaluate(const vertex<dimension>& v) const = 0;
};

}

#endif
