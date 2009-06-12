#ifndef SIMPLE_CFD_TRIANGULAR_CELL_HPP
#define SIMPLE_CFD_TRIANGULAR_CELL_HPP

#include<vector>
#include<cassert>
#include<map>
#include<boost/array.hpp>
#include<ostream>
#include<iostream>
#include<algorithm>
#include<numeric>
#include<cstddef>
#include"simple_cfd_fwd.hpp"
#include"vertex.hpp"
#include"mesh_topology.hpp"
#include"general_cell.hpp"

namespace cfd
{

class TriangularCell : public GeneralCell
{
public:
  static const std::size_t dimension = 2;
  typedef vertex<dimension> vertex_type;
  static const unsigned int vertex_count = 3;

  TriangularCell();
  virtual std::size_t getDimension() const;
  virtual std::size_t getVerticesPerCell() const;
  static std::map<vertex_type, double> getReferenceQuadrature();
  virtual std::map<vertex_type, double> getQuadrature(const mesh<TriangularCell>& m, const MeshEntity& entity) const;
  double getArea(const mesh<TriangularCell>& m, const MeshEntity& entity) const;
  vertex_type reference_to_physical(const mesh<TriangularCell>& m, const std::size_t cid, const vertex_type& vertex) const;
  bool contains(const mesh<TriangularCell>& m, const std::size_t cid, const vertex_type& v) const;
  std::set< std::set<std::size_t> > getIncidentVertices(MeshTopology& topology, const MeshEntity& cellEntity, std::size_t d) const;
};

}
#endif
