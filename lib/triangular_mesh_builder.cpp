#include <simple_cfd_fwd.hpp>
#include <triangular_mesh_builder.hpp>
#include <mesh.hpp>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

namespace cfd
{

TriangularMeshBuilder::TriangularMeshBuilder(const double _width, const double _height, const double _maxCellArea) : 
  width(_width), height(_height), maxCellArea(_maxCellArea)
{
  assert(width > 0.0);
  assert(height > 0.0);
}

mesh<TriangularMeshBuilder::cell_type> TriangularMeshBuilder::buildMesh() const
{
  const double sideLength = std::sqrt(2*maxCellArea);
  const int x_size = static_cast<int>(std::ceil(width / sideLength)) + 1;
  const int y_size = static_cast<int>(std::ceil(height / sideLength)) + 1;
  mesh<cell_type> m;

  assert(x_size > 1);
  assert(y_size > 1);

  // Create vertices
  vertex_id vid=0;

  const int x_nodes = x_size;
  const int y_nodes = y_size;

  for(int y=0; y < y_nodes; ++y)
  {
    for(int x=0; x < x_nodes; ++x)
    {
      m.addVertex(vid, vertex_type(static_cast<double>(x) / (x_nodes-1) * width, static_cast<double>(y) / (y_nodes-1) * height));
      ++vid;
    }
  }

  /*   Triangle Node Numbering

       LL        UR

       2           1     0
                    _____
       |\           \    |
       | \           \   |
       |  \           \  |  
       |   \           \ |
       |    \           \|
       ------           
       0     1           2
  */

  // Now create cells assuming first vertex has id 0
  std::vector<vertex_id> lower_left_vertices;
  lower_left_vertices.push_back(0);
  lower_left_vertices.push_back(1);
  lower_left_vertices.push_back(x_nodes);

  std::vector<vertex_id> upper_right_vertices;
  upper_right_vertices.push_back(1 + x_nodes);
  upper_right_vertices.push_back(x_nodes);
  upper_right_vertices.push_back(1);

  cell_id cid = 0;

  for(int y = 0; y < y_size - 1; ++y)
  {
    for(int x =0; x < x_size - 1; ++x)
    {
      std::vector<vertex_id> offset_lower_left_vertices(lower_left_vertices);
      std::vector<vertex_id> offset_upper_right_vertices(upper_right_vertices);

      const int offset = x_nodes * y + x;

      std::transform(offset_lower_left_vertices.begin(), 
        offset_lower_left_vertices.end(), offset_lower_left_vertices.begin(), boost::lambda::_1 + offset);

      std::transform(offset_upper_right_vertices.begin(), 
        offset_upper_right_vertices.end(), offset_upper_right_vertices.begin(), boost::lambda::_1 + offset);

      m.addCell(cid, cell_type(offset_lower_left_vertices));
      ++cid;
      m.addCell(cid, cell_type(offset_upper_right_vertices));
      ++cid;
    }
  }

  return m;
}

}
