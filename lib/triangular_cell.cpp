#include <triangular_cell.hpp>
#include <mesh.hpp>
#include <vector>

namespace cfd
{

TriangularCell::TriangularCell()
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

std::map<TriangularCell::vertex_type, double> TriangularCell::getReferenceQuadrature()
{
  /* Cubic Gaussian quadrature values from "Finite Elements: A Gentle Introduction" by Henwood and Bonet */
  /* These co-ordinates are defined on the reference triangle {(0,0), (0,1), (0,1)} */ 
  std::map<vertex_type, double> weightings;
  weightings[vertex_type(0.0, 0.0)] = 3.0/120.0;
  weightings[vertex_type(1.0, 0.0)] = 3.0/120.0;
  weightings[vertex_type(0.0, 1.0)] = 3.0/120.0;
  weightings[vertex_type(0.5, 0.0)] = 8.0/120.0;
  weightings[vertex_type(0.5, 0.5)] = 8.0/120.0;
  weightings[vertex_type(0.0, 0.5)] = 8.0/120.0;
  weightings[vertex_type(1.0/3, 1.0/3)] = 27.0/120.0;
  return weightings;
}

std::map<TriangularCell::vertex_type, double> TriangularCell::getQuadrature(const mesh<TriangularCell>& m, const MeshEntity& entity) const
{
  const std::map<vertex_type, double> referenceWeightings(getReferenceQuadrature());
  const double scaling = getArea(m, entity) / 0.5; // 0.5 is area of reference triangle
  std::map<vertex_type, double> weightings;

  for(std::map<vertex_type, double>::const_iterator refIter(referenceWeightings.begin()); refIter!=referenceWeightings.end(); ++refIter)
    weightings[reference_to_physical(m, entity.getIndex(), refIter->first)] = refIter->second * scaling;

  return weightings;
}

double TriangularCell::getArea(const mesh<TriangularCell>& m, const MeshEntity& entity) const
{
  const std::vector<vertex_type> vertices(m.getCoordinates(entity.getIndex()));
  const double doubleArea = vertices[0][0] * (vertices[1][1] - vertices[2][1]) +
                            vertices[1][0] * (vertices[2][1] - vertices[0][1]) +
                            vertices[2][0] * (vertices[0][1] - vertices[1][1]);
  return doubleArea / 2.0;
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

bool TriangularCell::contains(const mesh<TriangularCell>& m, const std::size_t cid, const vertex_type& v) const
{
  const std::vector<vertex_type> vertices(m.getCoordinates(cid));
  bool contained = true;

  for(unsigned vid=0; vid<vertex_count; ++vid)
  {
    const unsigned v1_ind = vid;
    const unsigned v2_ind = (vid+1)%vertex_count;
    const unsigned v3_ind = (vid+2)%vertex_count;

    std::vector<double> edgeVector;
    std::transform(vertices[v1_ind].begin(), vertices[v1_ind].end(), vertices[v2_ind].begin(), std::back_inserter(edgeVector), std::minus<double>());

    // We calculate the perpendicular to the edge which gives us the coefficients for the place
    std::vector<double> edgeNormal(edgeVector.size());
    edgeNormal[0] = -edgeVector[1];
    edgeNormal[1] = edgeVector[0];

    const double d = std::inner_product(edgeNormal.begin(), edgeNormal.end(), vertices[v1_ind].begin(), 0.0);
    const double vDotProd = std::inner_product(edgeNormal.begin(), edgeNormal.end(), v.begin(), 0.0) - d;
    const double v3DotProd = std::inner_product(edgeNormal.begin(), edgeNormal.end(), vertices[v3_ind].begin(), 0.0) - d;

    if ((vDotProd<0 && v3DotProd>0) || (vDotProd>0 && v3DotProd<0))
      contained = false;
  }
  return contained;
}

std::set< std::set<std::size_t> > TriangularCell::getIncidentVertices(MeshTopology& topology, const MeshEntity& cellEntity, std::size_t d) const
{
  //NOTE: we rely on the fact that the vertices are sorted, otherwise this method
  //      would be non-deterministic
  assert(cellEntity.getDimension() == dimension);
  assert(d <= dimension);

  std::vector<std::size_t> vertices(topology.getIndices(cellEntity, 0));
  assert(vertices.size() == vertex_count);

  std::set< std::set<std::size_t> > result;

  if (d == 2)
  {
    // All vertices are incident to the cell
    result.insert(std::set<std::size_t>(vertices.begin(), vertices.end()));
  }
  else if (d == 1)
  {
    for(unsigned edge=0; edge<3; ++edge)
    {
      std::set<std::size_t> edgeVertexSet;
      edgeVertexSet.insert(vertices[edge]);
      edgeVertexSet.insert(vertices[(edge+1)%3]);
      result.insert(edgeVertexSet);
    }
  }
  else if (d == 0)
  {
    // Vertices are only incident to themselves
    for(unsigned i=0; i<vertices.size(); ++i)
    {
      std::set<std::size_t> singleVertexSet;
      singleVertexSet.insert(vertices[i]);
      result.insert(singleVertexSet);
    }
  }

  return result;
}

}
