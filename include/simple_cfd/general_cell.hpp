#ifndef SIMPLE_CFD_GENERIC_CELL_HPP
#define SIMPLE_CFD_GENERIC_CELL_HPP

#include <set>
#include <map>
#include <vector>
#include <cstddef>
#include <boost/array.hpp>
#include "simple_cfd_fwd.hpp"
#include "mesh_cell.hpp"

namespace cfd
{

template<std::size_t D>
class GeneralCell : public MeshCell
{
public:
  static const std::size_t dimension = D;

  virtual std::size_t numEntities(std::size_t dimension) const = 0;
  virtual std::size_t getLocalIndex(MeshTopology& topology, const std::size_t cid, const MeshEntity& entity) const = 0;
  virtual vertex<dimension> getLocalVertex(const std::size_t index) const = 0;
  virtual QuadraturePoints<dimension> getQuadrature(const boost::array<std::size_t, dimension>& degrees) const = 0;
  virtual double getArea(const CellVertices<dimension>& vertices) const = 0;
  virtual double getJacobian(const CellVertices<dimension>& vertices, const MeshEntity& localEntity, const vertex<dimension>& v) const = 0;
  virtual vertex<dimension> referenceToPhysical(const CellVertices<dimension>& vertices, const vertex<dimension>& vertex) const = 0;
  virtual GlobalTransformation<dimension, dimension> getLocalGlobalTransformation() const = 0;
  virtual LocalTransformation<dimension, dimension> getCellReferenceLocalTransformation() const = 0;
  virtual LocalTransformation<dimension-1, dimension> getFacetReferenceLocalTransformation(const std::size_t fid) const = 0;
  virtual Tensor<dimension> getFacetNormal(const CellVertices<dimension>& vertices, const std::size_t fid,
    const vertex<dimension>& v) const = 0;
  virtual const FiniteElement<dimension>& getCoordinateMapping() const = 0;
  virtual ~GeneralCell() {}
};

}

#endif
