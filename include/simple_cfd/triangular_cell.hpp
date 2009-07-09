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

class TriangularCell : public GeneralCell
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
  virtual std::size_t getVerticesPerCell() const;
  virtual vertex<2> getLocalVertex(const std::size_t index) const;
  virtual std::size_t getLocalIndex(MeshTopology& topology, std::size_t cid, const MeshEntity& entity) const;
  virtual QuadraturePoints<2> getQuadrature(const std::size_t degree) const;
  double getArea(const mesh<TriangularCell>& m, const MeshEntity& entity) const;
  double getJacobian(const mesh<TriangularCell>& m, const MeshEntity& entity, const vertex_type& v) const;
  vertex_type reference_to_physical(const mesh<TriangularCell>& m, const std::size_t cid, const vertex_type& vertex) const;
  std::vector< std::set<std::size_t> > getIncidentVertices(MeshTopology& topology, const MeshEntity& cellEntity, std::size_t d) const;
  Tensor<dimension, 1, double> getFacetNormal(const mesh<TriangularCell>& m, const std::size_t cid, const std::size_t fid, const vertex_type& v) const;
};

}
#endif
