#include <simple_cfd_fwd.hpp>
#include <triangular_mesh_builder.hpp>
#include <mesh.hpp>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <boost/shared_array.hpp>
#include <boost/lexical_cast.hpp>
#include <libtriangle.hpp>

namespace cfd
{

namespace {

void clear(triangulateio& t)
{
  memset(&t, 0, sizeof(triangulateio));
}

boost::shared_array<char> stringToCharArray(const std::string& s)
{
  boost::shared_array<char> out(new char[s.size() + 1]);
  strcpy(&out[0], s.c_str());
  return out;
}

void trifreeMembers(triangulateio& t)
{
  trifree(t.pointlist);
  trifree(t.pointattributelist);
  trifree(t.pointmarkerlist);
  trifree(t.trianglelist);
  trifree(t.triangleattributelist);
  trifree(t.trianglearealist);
  trifree(t.neighborlist);
  trifree(t.segmentlist);
  trifree(t.segmentmarkerlist);
  trifree(t.holelist);
  trifree(t.regionlist);
  trifree(t.edgelist);
  trifree(t.edgemarkerlist);
  trifree(t.normlist);
}

}

TriangularMeshBuilder::TriangularMeshBuilder(const double _width, const double _height, const double _maxCellArea) : 
  width(_width), height(_height), maxCellArea(_maxCellArea)
{
  assert(width > 0.0);
  assert(height > 0.0);
}

mesh<TriangularMeshBuilder::cell_type> TriangularMeshBuilder::buildMesh() const
{
  return buildMeshTriangle();
}

mesh<TriangularMeshBuilder::cell_type> TriangularMeshBuilder::buildMeshTriangle() const
{
  triangulateio in, out;
  clear(in);
  clear(out);

  in.numberofpoints = 4;
  in.numberofpointattributes = 0;
  in.pointlist = new double[in.numberofpoints * 2];
  
  // BL
  in.pointlist[0] = 0.0;
  in.pointlist[1] = 0.0;
  // TL
  in.pointlist[2] = 0.0;
  in.pointlist[3] = height;
  // TR
  in.pointlist[4] = width;
  in.pointlist[5] = height;
  // BR
  in.pointlist[6] = width;
  in.pointlist[7] = 0.0;

  in.numberofsegments = 4;
  in.segmentlist = new int[in.numberofsegments * 2];
  in.segmentlist[0] = 0;
  in.segmentlist[1] = 1;
  in.segmentlist[2] = 1;
  in.segmentlist[3] = 2;
  in.segmentlist[4] = 2;
  in.segmentlist[5] = 3;
  in.segmentlist[6] = 3;
  in.segmentlist[7] = 0;

  // Q=quiet, z=zero-indexing, a=area-constraint
  boost::shared_array<char> options = stringToCharArray(std::string("Qza") + 
    boost::lexical_cast<std::string>(maxCellArea));

  triangulate(options.get(), &in, &out, NULL);

  delete[] in.pointlist;
  delete[] in.segmentlist;

  mesh<TriangularMeshBuilder::cell_type> m;

  for(vertex_id vid = 0; vid < static_cast<unsigned>(out.numberofpoints); ++vid)
  {
    const vertex_type v(out.pointlist[vid*2], out.pointlist[vid*2+1]);
    m.addVertex(v);
  }

  std::vector<vertex_id> cellVertices(3);

  for(cell_id cid = 0; cid < static_cast<unsigned>(out.numberoftriangles); ++cid)
  {
    cellVertices[0] = out.trianglelist[cid*3];
    cellVertices[1] = out.trianglelist[cid*3+1];
    cellVertices[2] = out.trianglelist[cid*3+2];
    m.addCell(cid, cell_type(cellVertices));
  }

  trifreeMembers(out);
  return m;
}

mesh<TriangularMeshBuilder::cell_type> TriangularMeshBuilder::buildMeshOld() const
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
      m.addVertex(vertex_type(static_cast<double>(x) / (x_nodes-1) * width, static_cast<double>(y) / (y_nodes-1) * height));
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
