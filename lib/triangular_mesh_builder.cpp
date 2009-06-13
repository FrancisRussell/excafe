#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <utility>
#include <iterator>
#include <simple_cfd_fwd.hpp>
#include <triangular_mesh_builder.hpp>
#include <mesh.hpp>
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

void TriangularMeshBuilder::addPolygon(const Polygon& polygon, const int label)
{
  polygons.push_back(std::make_pair(polygon, label));
}

mesh<TriangularMeshBuilder::cell_type> TriangularMeshBuilder::buildMeshTriangle() const
{
  triangulateio in, out;
  clear(in);
  clear(out);

  std::vector<double> pointList;
  std::vector<int> segmentList;
  std::vector<int> segmentMarkerList;

  // BL
  pointList.push_back(0.0);
  pointList.push_back(0.0);
  // BR
  pointList.push_back(width);
  pointList.push_back(0.0);
  // TR
  pointList.push_back(width);
  pointList.push_back(height);
  // TL
  pointList.push_back(0.0);
  pointList.push_back(height);

  segmentList.push_back(0);
  segmentList.push_back(1);
  segmentList.push_back(1);
  segmentList.push_back(2);
  segmentList.push_back(2);
  segmentList.push_back(3);
  segmentList.push_back(3);
  segmentList.push_back(0);

  segmentMarkerList.push_back(1);
  segmentMarkerList.push_back(2);
  segmentMarkerList.push_back(3);
  segmentMarkerList.push_back(4);

  handlePolygons(pointList, segmentList, segmentMarkerList);

  in.pointlist = &pointList[0];
  in.numberofpoints = pointList.size() / 2;

  in.segmentlist = &segmentList[0];
  in.numberofsegments = segmentList.size() / 2;

  assert(static_cast<int>(segmentMarkerList.size()) == in.numberofsegments);
  in.segmentmarkerlist = &segmentMarkerList[0];

  // Q=quiet, e=output-edges, p=read-segments, z=zero-indexing, a=area-constraint
  boost::shared_array<char> options = stringToCharArray(std::string("Qepza") + 
    boost::lexical_cast<std::string>(maxCellArea));

  triangulate(options.get(), &in, &out, NULL);

  mesh<TriangularMeshBuilder::cell_type> m;

  for(vertex_id vid = 0; vid < static_cast<unsigned>(out.numberofpoints); ++vid)
  {
    const vertex_type v(out.pointlist[vid*2], out.pointlist[vid*2+1]);
    const vertex_id givenVid = m.addVertex(v);
    assert(vid == givenVid);
  }

  std::vector<vertex_id> cellVertices(3);

  for(cell_id cid = 0; cid < static_cast<unsigned>(out.numberoftriangles); ++cid)
  {
    cellVertices[0] = out.trianglelist[cid*3];
    cellVertices[1] = out.trianglelist[cid*3+1];
    cellVertices[2] = out.trianglelist[cid*3+2];
    const cell_id givenCid = m.addCell(cellVertices);
    assert(givenCid == cid);
  }

  m.finish();

  // We can only start to map edge markers back to edges after the mesh is able to 
  // determine topological relations.

  const std::size_t facetDim = m.getDimension() - 1;
  MeshFunction<int> facetNumbering(facetDim);

  for(int edge=0; edge < out.numberofedges; ++edge)
  {
    const std::size_t v1 = out.edgelist[edge*2];
    const std::size_t v2 = out.edgelist[edge*2+1];
    const int label = out.edgemarkerlist[edge];

    const MeshEntity v1Entity(0, v1);
    bool foundEdge = false;

    for(mesh<cell_type>::local_iterator facetIter = m.local_begin(v1Entity, facetDim); facetIter!=m.local_end(v1Entity, facetDim); ++facetIter)
    {
      const std::vector<std::size_t> vertexIndices(m.getIndices(*facetIter, 0));
      if (std::find(vertexIndices.begin(), vertexIndices.end(), v2) != vertexIndices.end())
      {
        foundEdge = true;
        facetNumbering(*facetIter) = label;
      }
    }
    assert(foundEdge);
  }

  m.setFacetLabelling(facetNumbering);
  trifreeMembers(out);
  return m;
}

void TriangularMeshBuilder::handlePolygons(std::vector<double>& pointList, 
  std::vector<int>& segmentList, std::vector<int>& segmentMarkerList) const
{
  for(std::vector< std::pair<Polygon, int> >::const_iterator polyIter(polygons.begin());
    polyIter!=polygons.end(); ++polyIter)
  {
    const Polygon polygon(polyIter->first);
    const int label = polyIter->second;

    std::vector<double> points;
    std::vector<int> segments;

    for(std::size_t node=0; node<polygon.getNumSides(); ++node)
    {
      const double angle = polygon.getRotation() + (2*M_PI)*node/polygon.getNumSides();
      const double x = polygon.getOrigin()[0] + sin(angle) * polygon.getRadius();
      const double y = polygon.getOrigin()[1] + cos(angle) * polygon.getRadius();
      points.push_back(x);
      points.push_back(y);
    }

    for(std::size_t edge=0; edge<polygon.getNumSides(); ++edge)
    {
      segments.push_back(edge%polygon.getNumSides() + pointList.size()/2);
      segments.push_back((edge+1)%polygon.getNumSides() + pointList.size()/2);
    }

    pointList.insert(pointList.end(), points.begin(), points.end());
    segmentList.insert(segmentList.end(), segments.begin(), segments.end());
    std::fill_n(std::back_inserter(segmentMarkerList), polygon.getNumSides(), label);
  }
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
      const vertex_id givenVid = m.addVertex(vertex_type(static_cast<double>(x) / (x_nodes-1) * width, static_cast<double>(y) / (y_nodes-1) * height));
      assert(vid == givenVid);
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

      const cell_id llCid = m.addCell(offset_lower_left_vertices);
      assert(llCid == cid);
      ++cid;
      const cell_id urCid = m.addCell(offset_upper_right_vertices);
      assert(urCid == cid);
      ++cid;
    }
  }

  m.finish();
  return m;
}

}
