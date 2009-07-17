#ifndef SIMPLE_CFD_BOUNDARY_CONDITION2_HPP
#define SIMPLE_CFD_BOUNDARY_CONDITION2_HPP

#include "simple_cfd_fwd.hpp"
#include <cstddef>

namespace cfd
{

template<std::size_t D, std::size_t R>
class BoundaryCondition2
{
public:
  static const std::size_t dimension = D;
  static const std::size_t rank = R;

  virtual bool applies(MeshTopology& topology, const MeshEntity& entity, const std::size_t label) const = 0;
  virtual bool constrainsIndex(MeshTopology& topology, const MeshEntity& entity, const std::size_t label) const = 0;
  virtual std::size_t constrainedIndex(MeshTopology& topology, const MeshEntity& entity, const std::size_t label) const = 0;
  virtual Tensor<dimension, rank> value(MeshTopology& topology, 
    const MeshEntity& entity, const std::size_t label, const vertex<dimension>& v) const = 0;
  virtual ~BoundaryCondition2() {}
};

}

#endif
