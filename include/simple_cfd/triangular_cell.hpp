#ifndef SIMPLE_CFD_TRIANGULAR_CELL_HPP
#define SIMPLE_CFD_TRIANGULAR_CELL_HPP

#include <vector>
#include <map>
#include <set>
#include <cstddef>
#include "simple_cfd_fwd.hpp"
#include "vertex.hpp"
#include "mesh_topology.hpp"
#include "general_cell.hpp"
#include "quadrature_points.hpp"

namespace cfd
{

class TriangularCell : public GeneralCell<2>
{
public:
  static const std::size_t dimension = 2;
  typedef vertex<dimension> vertex_type;
  static const unsigned int vertex_count = 3;

private:
  const std::vector<vertex_type> localVertices;

  static std::vector<vertex_type> buildLocalVertices();
  static std::map<vertex_type, double> buildCellQuadrature(const std::size_t degree);
  static std::map<vertex_type, double> normaliseQuadrature(const std::map<vertex_type, double>& quadrature, const double value); 

public:
  TriangularCell();
  virtual std::size_t getDimension() const;
  virtual std::size_t numEntities(std::size_t dimension) const;
  virtual vertex<dimension> getLocalVertex(const std::size_t index) const;
  virtual std::size_t getLocalIndex(MeshTopology& topology, std::size_t cid, const MeshEntity& entity) const;
  virtual QuadraturePoints<dimension> getQuadrature(const std::size_t degree) const;
  virtual double getArea(const CellVertices<dimension>& vertices) const;
  double getJacobian(const CellVertices<dimension>& vertices, const MeshEntity& localEntity, const vertex_type& v) const;
  vertex_type referenceToPhysical(const CellVertices<dimension>& vertices, const vertex_type& vertex) const;
  std::vector< std::set<std::size_t> > getIncidentVertices(MeshTopology& topology, const MeshEntity& cellEntity, std::size_t d) const;
  Tensor<dimension, 1, double> getFacetNormal(const CellVertices<dimension>& vertices, const std::size_t fid, const vertex_type& v) const;
  virtual std::auto_ptr<MeshCell> cloneMeshCell() const;
  virtual std::auto_ptr< GeneralCell<dimension> > cloneGeneralCell() const;
};

}
#endif
