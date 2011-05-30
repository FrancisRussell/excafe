#ifndef SIMPLE_CFD_TRIANGULAR_CELL_HPP
#define SIMPLE_CFD_TRIANGULAR_CELL_HPP

#include <vector>
#include <map>
#include <set>
#include <cstddef>
#include <boost/array.hpp>
#include <boost/scoped_ptr.hpp>
#include "simple_cfd_fwd.hpp"
#include "vertex.hpp"
#include "mesh_topology.hpp"
#include "general_cell.hpp"
#include "quadrature_points.hpp"
#include "finite_element.hpp"
#include "coordinate_transformation.hpp"

namespace cfd
{

class TriangularCell : public GeneralCell<2>
{
public:
  static const std::size_t dimension = 2;
  typedef vertex<dimension> vertex_type;
  static const unsigned int vertex_count = 3;

private:
  friend class util::Singleton<TriangularCell>;
  const std::vector<vertex_type> localVertices;
  mutable boost::scoped_ptr< FiniteElement<dimension> > coordinateMapping;

  static std::vector<vertex_type> buildLocalVertices();
  static std::map<vertex_type, double> buildCellQuadrature(const boost::array<std::size_t, dimension>& degrees);
  static std::map<vertex_type, double> normaliseQuadrature(const std::map<vertex_type, double>& quadrature, const double value); 

  TriangularCell();

public:
  virtual std::size_t getDimension() const;
  virtual std::size_t numEntities(std::size_t dimension) const;
  virtual vertex<dimension> getLocalVertex(const std::size_t index) const;
  virtual std::size_t getLocalIndex(MeshTopology& topology, std::size_t cid, const MeshEntity& entity) const;
  virtual QuadraturePoints<dimension> getQuadrature(const boost::array<std::size_t, dimension>& degrees) const;
  virtual double getArea(const CellVertices<dimension>& vertices) const;
  double getJacobian(const CellVertices<dimension>& vertices, const MeshEntity& localEntity, const vertex_type& v) const;
  vertex_type referenceToPhysical(const CellVertices<dimension>& vertices, const vertex_type& vertex) const;
  virtual GlobalTransformation<dimension, dimension> getLocalGlobalTransformation() const;
  virtual LocalTransformation<dimension, dimension> getCellReferenceLocalTransformation() const;
  virtual LocalTransformation<dimension-1, dimension> getFacetReferenceLocalTransformation(const std::size_t fid) const;
  virtual std::set<std::size_t> getIncidentVertices(const MeshEntity& localEntity) const;
  virtual std::set<std::size_t> getIncidentVertices(MeshTopology& topology, const std::size_t cid, const MeshEntity& localEntity) const;
  Tensor<dimension> getFacetNormal(const CellVertices<dimension>& vertices, const std::size_t fid, const vertex_type& v) const;
  virtual const FiniteElement<dimension>& getCoordinateMapping() const;
};

}
#endif
