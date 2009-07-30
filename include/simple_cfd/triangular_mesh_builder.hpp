#ifndef SIMPLE_CFD_TRIANGULAR_MESH_BUILDER_HPP
#define SIMPLE_CFD_TRIANGULAR_MESH_BUILDER_HPP

#include <simple_cfd_fwd.hpp>
#include <mesh.hpp>
#include <vector>
#include <utility>
#include <libtriangle.hpp>
#include <polygon.hpp>

namespace cfd
{

class TriangularMeshBuilder
{
private:
  typedef TriangularCell cell_type;
  typedef TriangularCell::vertex_type vertex_type;
  static const std::size_t dimension = cell_type::dimension;
  const double width;
  const double height;
  const double maxCellArea;
  std::vector< std::pair<Polygon, int> > polygons;

  Mesh<dimension> buildMeshOld() const;
  Mesh<dimension> buildMeshTriangle() const;
  void handlePolygons(std::vector<double>& pointList, 
    std::vector<int>& segmentList, std::vector<int>& segmentMarkerList,
    std::vector<double>& holeList) const;

public:
  TriangularMeshBuilder(const double width, const double height, const double maxCellArea);
  void addPolygon(const Polygon& polygon, int label);
  Mesh<dimension> buildMesh() const;
};

}
#endif
