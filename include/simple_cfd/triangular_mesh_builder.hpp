#ifndef SIMPLE_CFD_TRIANGULAR_MESH_BUILDER_HPP
#define SIMPLE_CFD_TRIANGULAR_MESH_BUILDER_HPP

#include <simple_cfd_fwd.hpp>
#include <mesh.hpp>
#include <libtriangle.hpp>

namespace cfd
{

class TriangularMeshBuilder
{
private:
  typedef TriangularCell cell_type;
  typedef TriangularCell::vertex_type vertex_type;
  const double width;
  const double height;
  const double maxCellArea;

  mesh<cell_type> buildMeshOld() const;
  mesh<cell_type> buildMeshTriangle() const;

public:
  TriangularMeshBuilder(const double width, const double height, const double maxCellArea);
  mesh<TriangularCell> buildMesh() const;
};

}
#endif
