#ifndef EXCAFE_FUNCTION_HPP
#define EXCAFE_FUNCTION_HPP

#include "excafe_fwd.hpp"

namespace excafe
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
