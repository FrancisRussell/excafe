#include <vector>
#include <cassert>
#include <map>
#include <set>
#include <numeric>
#include <algorithm>
#include <utility>
#include <cmath>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <simple_cfd_fwd.hpp>
#include <triangular_cell.hpp>
#include <lagrange_triangle_linear.hpp>
#include <numeric/tensor.hpp>
#include <numeric/quadrature.hpp>
#include <symbolic/rational.hpp>
#include <mesh.hpp>
#include <exception.hpp>
#include <cell_vertices.hpp>

namespace cfd
{

TriangularCell::TriangularCell() : localVertices(buildLocalVertices())
{
}

std::size_t TriangularCell::getDimension() const
{
  return dimension;
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

std::map<TriangularCell::vertex_type, double> TriangularCell::buildCellQuadrature(const boost::array<std::size_t, dimension>& degrees)
{
  Quadrature quadrature;

  // We increase degree of s quadrature by one to handle (1-s) factor in Jacobian of co-ordinate
  // transformation from reference square to triangle.
  const std::map<double, double> rQuadrature(quadrature.getGauss(degrees[0]));
  const std::map<double, double> sQuadrature(quadrature.getGauss(degrees[1]+1));

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

QuadraturePoints<2> TriangularCell::getQuadrature(const boost::array<std::size_t, dimension>& degrees) const
{
  Quadrature quadrature;
  const std::map<double, double> unitQuadrature = quadrature.getGauss(std::max(degrees[0], degrees[1]));

  QuadraturePoints<2> points;
  points.setQuadrature(MeshEntity(dimension, 0), buildCellQuadrature(degrees));

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

double TriangularCell::getArea(const CellVertices<2>& vertices) const
{
  // Jacobian is constant so v doesn't matter
  const vertex_type v(0.0, 0.0);
  const double area = 0.5 * getJacobian(vertices, MeshEntity(dimension, 0), v);
  assert(area >= 0.0);
  return area;
}

double TriangularCell::getJacobian(const CellVertices<2>& vertices, const MeshEntity& localEntity, const vertex_type& b) const
{
  if (localEntity.getDimension() == 2)
  {
    assert(vertices.size() == 3);
    const double jacobian = vertices[0][0] * (vertices[1][1] - vertices[2][1]) +
                            vertices[1][0] * (vertices[2][1] - vertices[0][1]) +
                            vertices[2][0] * (vertices[0][1] - vertices[1][1]);
    return jacobian;
  }
  else if (localEntity.getDimension() == 1)
  {
    const std::size_t localIndex = localEntity.getIndex();
    const vertex_type difference = vertices[localIndex] - vertices[(localIndex+1)%3];
    return std::sqrt(difference[0] * difference[0] + difference[1] * difference[1]);
  }
  else
  {
    assert(false);
  }

  return 0.0;
}

TriangularCell::vertex_type TriangularCell::referenceToPhysical(const CellVertices<dimension>& vertices, const vertex_type& vertex) const
{
  const double xsi = vertex[0];
  const double eta = vertex[1];

  const vertex_type v = (1.0 - xsi - eta) * vertices[0] +
                        xsi * vertices[1] +
                        eta *vertices[2];

  return v;
}

std::set<std::size_t> TriangularCell::getIncidentVertices(const MeshEntity& localEntity) const
{
  const std::size_t d = localEntity.getDimension();
  const std::size_t i = localEntity.getIndex();

  assert(d <= dimension);
  assert(i < numEntities(d));

  std::set<std::size_t> result;

  if (d == 2)
  {
    // All vertices are incident to the cell
    for(std::size_t v=0; v<vertex_count; ++v)
      result.insert(v);
  }
  else if (d == 1)
  {
    result.insert((i+0)%3);
    result.insert((i+1)%3);
  }
  else if (d == 0)
  {
    // Vertices are only incident to themselves
    result.insert(i);
  }

  return result;
}

std::set<std::size_t> TriangularCell::getIncidentVertices(MeshTopology& topology, const std::size_t cid, const MeshEntity& localEntity) const
{
  const MeshEntity cellEntity(dimension, cid);
  const std::vector<std::size_t> vertices(topology.getIndices(cellEntity, 0));
  assert(vertices.size() == vertex_count);

  const std::set<std::size_t> localIndices = getIncidentVertices(localEntity);
  std::set<std::size_t> globalIndices;

  BOOST_FOREACH(const std::size_t localIndex, localIndices)
    globalIndices.insert(vertices[localIndex]);

  return globalIndices;
}

std::size_t TriangularCell::getLocalIndex(MeshTopology& topology, const std::size_t cid, const MeshEntity& entity) const
{
  // Get vertices on entity
  const std::vector<std::size_t> vertexIndices(topology.getIndices(entity, 0));
  const std::set<std::size_t> vertexIndicesSet(vertexIndices.begin(), vertexIndices.end());

  for(std::size_t e=0; e<numEntities(entity.getDimension()); ++e)
  {
    const MeshEntity localEntity(entity.getDimension(), e);
    const std::set<std::size_t> incident(getIncidentVertices(topology, cid, localEntity));

    if (incident == vertexIndicesSet)
      return e;
  }

  CFD_EXCEPTION("Requested entity not present on cell.");
}

Tensor<TriangularCell::dimension> TriangularCell::getFacetNormal(const CellVertices<2>& vertices, const std::size_t localFacetID, const vertex_type& v) const
{
  Tensor<dimension> normal(1);
  assert(vertices.size() == 3);

  const vertex_type v1 = vertices[localFacetID];
  const vertex_type v2 = vertices[(localFacetID+1)%3];

  const double facetLength = v1.distance(v2);
  const vertex_type delta = v2 - v1;

  // The x and y values are exchanged to create the perpendicular
  normal(0) = delta[1] / facetLength;
  normal(1) = -delta[0] / facetLength;

  return normal;
}

TriangularCell::vertex_type TriangularCell::getLocalVertex(const std::size_t index) const
{
  std::vector<vertex_type::value_type> components;

  BOOST_FOREACH(const exact_vertex_type::value_type& c, getLocalVertexExact(index))
    components.push_back(c.toFloat().toDouble());

  return vertex_type(components);
}

TriangularCell::exact_vertex_type TriangularCell::getLocalVertexExact(const std::size_t index) const
{
  assert(index < localVertices.size());
  return localVertices[index];
}

std::vector<TriangularCell::exact_vertex_type> TriangularCell::buildLocalVertices()
{
  std::vector<exact_vertex_type> vertices;
  vertices.push_back(exact_vertex_type(0, 0));
  vertices.push_back(exact_vertex_type(1, 0));
  vertices.push_back(exact_vertex_type(0, 1));
  return vertices;
}

std::size_t TriangularCell::numEntities(const std::size_t d) const
{
  assert(d <= dimension);
  std::size_t numEntitiesArray[dimension+1] = {3, 3, 1};
  return numEntitiesArray[d];
}

const FiniteElement<TriangularCell::dimension>& TriangularCell::getCoordinateMapping() const
{
  //IMPORTANT: we need this to be lazily constructed to avoid an initialisation dependence loop

  if (!coordinateMapping)
    coordinateMapping.reset(new LagrangeTriangleLinear<0>());

  return *coordinateMapping;
}

GlobalTransformation<2, 2> TriangularCell::getLocalGlobalTransformation() const
{
  return GlobalTransformation<dimension, dimension>(getCoordinateMapping(), *this);
}

LocalTransformation<2, 2> TriangularCell::getCellReferenceLocalTransformation() const
{
  detail::PositionPlaceholder position;
  SmallVector<2, LocalTransformation<1, 2>::expression_t> transform;

  transform[0] = (position[0]*0.5 + 0.5)*(0.5 - position[1]*0.5);
  transform[1] = position[1]*0.5 + 0.5;

  return LocalTransformation<2,2>(transform);
}

LocalTransformation<1, 2> TriangularCell::getFacetReferenceLocalTransformation(const std::size_t fid) const
{
  typedef LocalTransformation<1, 2>::expression_t expression_t;

  detail::PositionPlaceholder refPosition;
  detail::PositionPlaceholder localPosition;
  SmallVector<2, expression_t> transform(getCellReferenceLocalTransformation().getTransformed());
  expression_t::value_map valueMap;

  if (fid == 0)
  {
    valueMap.bind(localPosition[0], refPosition[0]);
    valueMap.bind(localPosition[1], -1.0);
  }
  else if (fid == 1)
  {
    valueMap.bind(localPosition[0], 1.0);
    valueMap.bind(localPosition[1], refPosition[0]);
  }
  else if (fid == 2)
  {
    valueMap.bind(localPosition[0], -1.0);
    valueMap.bind(localPosition[1], refPosition[0]);
  }
  else
  {
    CFD_EXCEPTION("Transform requested for invalid facet.");
  }

  for(std::size_t d=0; d<transform.numRows(); ++d)
  {
    transform[d].substituteValues(valueMap);
  }

  return transform;
}

std::vector<TriangularCell::exact_vertex_type>
TriangularCell::constructLattice(const exact_vertex_type& base, const std::vector<exact_vertex_type>& offsets, const std::size_t degree)
{
  // Note: this code deliberately only constructs interior points.

  std::vector<exact_vertex_type> result;

  if (offsets.empty())
  {
    result.push_back(base);
  }
  else
  {
    const exact_vertex_type offset = offsets.front();
    const std::vector<exact_vertex_type> subOffsets(offsets.begin()+1, offsets.end());

    for(std::size_t i=1; i<degree; ++i)
    {
      const exact_vertex_type newBase = base + offset*i;
      const std::vector<exact_vertex_type> subLattice = constructLattice(newBase, subOffsets, degree-i);
      result.insert(result.end(), subLattice.begin(), subLattice.end());
    }
  }
  
  return result;
}

std::vector<TriangularCell::exact_vertex_type> 
TriangularCell::getSimplexPoints(const std::vector<exact_vertex_type>& vertices, const std::size_t degree)
{
  assert(!vertices.empty());

  const exact_vertex_type base = vertices.front();
  std::vector<exact_vertex_type> offsets;

  BOOST_FOREACH(const exact_vertex_type& v, std::make_pair(vertices.begin() + 1, vertices.end()))
    offsets.push_back((v - base)/degree);

  return constructLattice(base, offsets, degree);
}

std::map< MeshEntity, std::vector<TriangularCell::exact_vertex_type> > 
TriangularCell::getPoints(const std::size_t degree) const
{
  if (degree < 1)
    CFD_EXCEPTION("Lattice degree must be greater than 0.");

  std::map< MeshEntity, std::vector<exact_vertex_type> > pointMap;

  for(std::size_t d=0; d<=dimension; ++d)
  {
    for(std::size_t e=0; e<numEntities(d); ++e)
    {
      const MeshEntity entity(d, e);
      const std::set<std::size_t> vertexIndices = getIncidentVertices(entity);

      std::vector<exact_vertex_type> vertices;
      BOOST_FOREACH(const std::size_t index, vertexIndices)
        vertices.push_back(getLocalVertexExact(index));

      pointMap[entity] = getSimplexPoints(vertices, degree);
    }
  }

  return pointMap;
}

}
