#ifndef EXCAFE_GENERIC_CELL_HPP
#define EXCAFE_GENERIC_CELL_HPP

#include <set>
#include <map>
#include <vector>
#include <cstddef>
#include <boost/array.hpp>
#include "excafe_fwd.hpp"
#include "symbolic/rational.hpp"
#include "mesh_cell.hpp"

namespace excafe
{

template<std::size_t D>
class GeneralCell : public MeshCell
{
public:
  static const std::size_t dimension = D;
  typedef vertex<dimension> vertex_type;
  typedef vertex<dimension, symbolic::Rational> exact_vertex_type;

  virtual std::size_t getLocalIndex(MeshTopology& topology, const std::size_t cid, const MeshEntity& entity) const = 0;
  virtual vertex_type getLocalVertex(const std::size_t index) const = 0;
  virtual QuadraturePoints<dimension> getQuadrature(const boost::array<std::size_t, dimension>& degrees) const = 0;
  virtual double getArea(const CellVertices<dimension>& vertices) const = 0;
  virtual double getJacobian(const CellVertices<dimension>& vertices, const MeshEntity& localEntity, const vertex_type& v) const = 0;
  virtual vertex_type referenceToPhysical(const CellVertices<dimension>& vertices, const vertex_type& vertex) const = 0;
  virtual GlobalTransformation<dimension, dimension> getLocalGlobalTransformation() const = 0;
  virtual LocalTransformation<dimension, dimension> getCellReferenceLocalTransformation() const = 0;
  virtual LocalTransformation<dimension-1, dimension> getFacetReferenceLocalTransformation(const std::size_t fid) const = 0;
  virtual Tensor<dimension> getFacetNormal(const CellVertices<dimension>& vertices, const std::size_t fid,
    const vertex_type& v) const = 0;
  virtual const FiniteElement<dimension>& getCoordinateMapping() const = 0;
  virtual std::map< MeshEntity, std::vector<exact_vertex_type> > getPoints(const std::size_t degree) const = 0;
  virtual ~GeneralCell() {}
};

}

#endif
