#include <vector>
#include <cassert>
#include <map>
#include <set>
#include <numeric>
#include <algorithm>
#include <utility>
#include <cmath>
#include <simple_cfd_fwd.hpp>
#include <triangular_cell.hpp>
#include <numeric/tensor.hpp>
#include <numeric/quadrature.hpp>
#include <mesh.hpp>

namespace cfd
{

TriangularCell::TriangularCell() : localVertices(buildLocalVertices())
{
}

std::size_t TriangularCell::getDimension() const
{
  return dimension;
}

std::size_t TriangularCell::getVerticesPerCell() const
{
  return vertex_count;
}

std::map<TriangularCell::vertex_type, double> TriangularCell::normaliseQuadrature(const std::map<vertex_type, double>& quadrature, const double value)
{
  double sum = 0.0;
  std::map<vertex_type, double> newQuadrature;

  for(std::map<vertex_type, double>::const_iterator qIter(quadrature.begin()); qIter!=quadrature.end(); ++qIter)
    sum += qIter->second;

  for(std::map<vertex_type, double>::const_iterator qIter(quadrature.begin()); qIter!=quadrature.end(); ++qIter)
    newQuadrature.insert(std::make_pair(qIter->first, qIter->second * (value/sum)));

  return newQuadrature;
}

std::map<TriangularCell::vertex_type, double> TriangularCell::buildCellQuadrature(const std::size_t degree)
{
  Quadrature quadrature;

  // We increase degree of s quadrature by one to handle (1-s) factor in Jacobian of co-ordinate
  // transformation from reference square to triangle.
  const std::map<double, double> rQuadrature(quadrature.getGauss(degree));
  const std::map<double, double> sQuadrature(quadrature.getGauss(degree+1));

  std::map<vertex_type, double> squareQuadrature;

  for(std::map<double, double>::const_iterator rQuadIter(rQuadrature.begin()); rQuadIter!=rQuadrature.end(); ++rQuadIter)
    for(std::map<double, double>::const_iterator sQuadIter(sQuadrature.begin()); sQuadIter!=sQuadrature.end(); ++sQuadIter)
      squareQuadrature[vertex_type(rQuadIter->first, sQuadIter->first)] = rQuadIter->second * sQuadIter->second;

  std::map<vertex_type, double> triangularQuadrature;

  for(std::map<vertex_type, double>::const_iterator sIter(squareQuadrature.begin()); sIter!=squareQuadrature.end(); ++sIter)
  {
    const vertex_type oldLocation(sIter->first);
    const vertex_type newLocation((oldLocation[0]*0.5 + 0.5)*(0.5 - oldLocation[1]*0.5), oldLocation[1]*0.5 + 0.5);

    triangularQuadrature[newLocation] += sIter->second * (1-oldLocation[1]) / 8.0;
  }

  return triangularQuadrature;
}

QuadraturePoints<2> TriangularCell::getQuadrature(const std::size_t degree) const
{
  Quadrature quadrature;
  const std::map<double, double> unitQuadrature = quadrature.getGauss(degree);

  QuadraturePoints<2> points;
  points.setQuadrature(MeshEntity(dimension, 0), buildCellQuadrature(degree));

  for(std::size_t index=0; index<3; ++index)
  {
    std::map<vertex_type, double> facetWeightings;

    for(std::map<double, double>::const_iterator uIter(unitQuadrature.begin()); uIter!=unitQuadrature.end(); ++uIter)
    {
      double x = 0.0;
      double y = 0.0;

      //NOTE: the edge mapping MUST match reference_to_physical and getLocalIndex
      if (index == 0)
      {
        x = uIter->first;
        y = -1.0;
      }
      else if (index == 1)
      {
        x = 1.0;
        y = uIter->first;
      }
      else if (index == 2)
      {
        x = -1.0;
        y = uIter->first;
      }
      else
      {
        assert(false);
      }

      const vertex_type location((x*0.5 + 0.5)*(0.5 - y*0.5), y*0.5 + 0.5);
      facetWeightings[location] = uIter->second;
    }
    
    points.setQuadrature(MeshEntity(dimension-1, index), normaliseQuadrature(facetWeightings, 1.0));
  }

  return points;
}

double TriangularCell::getArea(const mesh<TriangularCell>& m, const MeshEntity& entity) const
{
  // Jacobian is constant so v doesn't matter
  const vertex_type v(0.0, 0.0);
  const double area = 0.5 * getJacobian(m, entity, v);
  assert(area >= 0.0);
  return area;
}

double TriangularCell::getJacobian(const mesh<TriangularCell>& m, const MeshEntity& entity, const vertex_type& b) const
{
  const std::vector<vertex_type> vertices(m.getCoordinates(m.getContainingCell(entity)));

  if (entity.getDimension() == 2)
  {
    assert(vertices.size() == 3);
    const double jacobian = vertices[0][0] * (vertices[1][1] - vertices[2][1]) +
                            vertices[1][0] * (vertices[2][1] - vertices[0][1]) +
                            vertices[2][0] * (vertices[0][1] - vertices[1][1]);
    return jacobian;
  }
  else if (entity.getDimension() == 1)
  {
    const std::size_t localIndex = getLocalIndex(m.getTopology(), m.getContainingCell(entity), entity);
    const vertex_type difference = vertices[localIndex] - vertices[(localIndex+1)%3];
    return std::sqrt(difference[0] * difference[0] + difference[1] * difference[1]);
  }
  else
  {
    assert(false);
  }

  return 0.0;
}

TriangularCell::vertex_type TriangularCell::reference_to_physical(const mesh<TriangularCell>& m, const std::size_t cid, const vertex_type& vertex) const
{
  const std::vector<vertex_type> vertices(m.getCoordinates(cid));

  const double xsi = vertex[0];
  const double eta = vertex[1];

  const vertex_type v = (1.0 - xsi - eta) * vertices[0] +
                        xsi * vertices[1] +
                        eta *vertices[2];

  return v;
}

std::vector< std::set<std::size_t> > TriangularCell::getIncidentVertices(MeshTopology& topology, const MeshEntity& cellEntity, std::size_t d) const
{
  assert(cellEntity.getDimension() == dimension);
  assert(d <= dimension);

  std::vector<std::size_t> vertices(topology.getIndices(cellEntity, 0));
  assert(vertices.size() == vertex_count);

  std::vector< std::set<std::size_t> > result;

  if (d == 2)
  {
    // All vertices are incident to the cell
    result.push_back(std::set<std::size_t>(vertices.begin(), vertices.end()));
  }
  else if (d == 1)
  {
    for(unsigned edge=0; edge<3; ++edge)
    {
      std::set<std::size_t> edgeVertexSet;
      edgeVertexSet.insert(vertices[edge]);
      edgeVertexSet.insert(vertices[(edge+1)%3]);
      result.push_back(edgeVertexSet);
    }
  }
  else if (d == 0)
  {
    // Vertices are only incident to themselves
    for(unsigned i=0; i<vertices.size(); ++i)
    {
      std::set<std::size_t> singleVertexSet;
      singleVertexSet.insert(vertices[i]);
      result.push_back(singleVertexSet);
    }
  }

  return result;
}

std::size_t TriangularCell::getLocalIndex(MeshTopology& topology, const std::size_t cid, const MeshEntity& entity) const
{
  // Get vertices on entity
  const std::vector<std::size_t> vertexIndices(topology.getIndices(entity, 0));
  const std::set<std::size_t> vertexIndicesSet(vertexIndices.begin(), vertexIndices.end());

  // Get sets of vertices corresponding to entities of dimension entity.getDimension() on cell cid.
  const std::vector< std::set<std::size_t> > incidentVertices(getIncidentVertices(topology, MeshEntity(dimension, cid), entity.getDimension()));
  
  const std::vector< std::set<std::size_t> >::const_iterator 
    entityIter(std::find(incidentVertices.begin(), incidentVertices.end(), vertexIndicesSet));

  // If this assertion fails, entity wasn't present on cell cid
  assert(entityIter != incidentVertices.end());

  return entityIter - incidentVertices.begin();
}


Tensor<TriangularCell::dimension, 1, double> TriangularCell::getFacetNormal(const mesh<TriangularCell>& m, 
  const std::size_t cid, const std::size_t fid, const vertex_type& v) const
{
  Tensor<dimension, 1, double> normal;
  std::vector< vertex<dimension> > cellVertices(m.getCoordinates(cid));
  assert(cellVertices.size() == 3);

  const std::size_t localFacetID = getLocalIndex(m.getTopology(), cid, MeshEntity(dimension-1, fid));
  const vertex_type v1 = cellVertices[localFacetID];
  const vertex_type v2 = cellVertices[(localFacetID+1)%3];

  const std::size_t zero = 0;
  const std::size_t one = 1;

  const double facetLength = v1.distance(v2);
  const vertex_type delta = v2 - v1;

  // The x and y values are exchanged to create the perpendicular
  normal[&zero] = delta[1] / facetLength;
  normal[&one] = -delta[0] / facetLength;

  return normal;
}

TriangularCell::vertex_type TriangularCell::getLocalVertex(const std::size_t index) const
{
  assert(index < localVertices.size());
  return localVertices[index];
}

std::vector<TriangularCell::vertex_type> TriangularCell::buildLocalVertices()
{
  std::vector<vertex_type> vertices;
  vertices.push_back(vertex_type(0.0, 0.0));
  vertices.push_back(vertex_type(1.0, 0.0));
  vertices.push_back(vertex_type(0.0, 1.0));
  return vertices;
}

}
