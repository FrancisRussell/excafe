#ifndef SIMPLE_CFD_FINITE_ELEMENT_HPP
#define SIMPLE_CFD_FINITE_ELEMENT_HPP

#include <vector>
#include <set>
#include <cstddef>
#include <map>
#include <boost/tuple/tuple.hpp>
#include "simple_cfd_fwd.hpp"
#include "dof_association.hpp"

namespace cfd
{

template<std::size_t D>
class FiniteElement
{
public:
  static const std::size_t dimension = D;
  typedef vertex<dimension> vertex_type;

  virtual unsigned getTensorIndex(const Mesh<dimension>& mesh, const std::size_t cid, const std::size_t dof) const = 0;
  virtual unsigned spaceDimension() const = 0; // Number of basis functions
  virtual std::vector< std::pair<unsigned, unsigned> > getCommonDegreesOfFreedom(const Mesh<dimension>& m, const cell_id cid, const cell_id cid2) const = 0;
  virtual vertex_type getDofCoordinateLocal(const std::size_t dof) const = 0;
  virtual vertex_type getDofCoordinateGlobal(const Mesh<dimension>& m, const cell_id cid, const std::size_t dof) const = 0;
  virtual std::set< boost::tuple<const FiniteElement<dimension>*, cell_id, std::size_t> > getDegreesOfFreedom(MeshTopology& topology, 
    const cell_id cid, const MeshEntity& entity) const = 0;
  virtual MeshEntity getLocalDofMeshAssociation(const Mesh<dimension>& mesh, const std::size_t cid, const std::size_t dof) const = 0;
  virtual ~FiniteElement() {}
};

}

#endif
