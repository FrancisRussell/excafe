#ifndef SIMPLE_CFD_BOUNDARY_CONDITION3_HPP
#define SIMPLE_CFD_BOUNDARY_CONDITION3_HPP

#include <cstddef>
#include "simple_cfd_fwd.hpp"

namespace cfd
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
