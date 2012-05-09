#ifndef EXCAFE_BOUNDARY_CONDITION3_HPP
#define EXCAFE_BOUNDARY_CONDITION3_HPP

#include <cstddef>
#include "excafe_fwd.hpp"

namespace excafe
{

template<std::size_t D>
class BoundaryCondition3
{
public:
  static const std::size_t dimension = D;

  virtual std::size_t getRank() const = 0;
  virtual bool applies(const int label) const = 0; 
  virtual Tensor<dimension> getValue(const vertex<dimension>& location, const int label) const = 0; 
  virtual ~BoundaryCondition3() {}
};

}

#endif
