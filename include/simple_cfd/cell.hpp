#ifndef SIMPLE_CFD_CELL_HPP
#define SIMPLE_CFD_CELL_HPP

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
#include"mesh.hpp"
#include"vertex.hpp"
#include"mesh_topology.hpp"
#include"general_cell.hpp"

namespace cfd
{

template<shape s>
class cell 
{
};

template<>
class cell<triangle> : public GeneralCell
{
public:
  static const shape cell_shape = triangle;
  static const unsigned int dimension = shape_dimensions<triangle>::dimension;
  typedef vertex<dimension> vertex_type;
  static const unsigned int vertex_count = 3;

private:
  boost::array<vertex_id, vertex_count> vertex_ids;

public:
  cell()
  {
  }

  cell(const std::vector<vertex_id>& vertices)
  {
    assert(vertices.size() == vertex_count);
    std::copy(vertices.begin(), vertices.end(), vertex_ids.begin());
  }

  virtual std::size_t getDimension() const
  {
    return dimension;
  }

  std::vector<vertex_id> getIndices() const
  {
    const std::vector<vertex_id> indices(vertex_ids.begin(), vertex_ids.end());
    return indices;
  }

  std::vector<vertex_type> getCoordinates(const mesh_geometry<dimension>& geometry) const
  {
    std::vector< vertex<dimension> > coords;

    for(boost::array<vertex_id, vertex_count>::const_iterator vertexIter(vertex_ids.begin()); vertexIter!=vertex_ids.end(); ++vertexIter)
      coords.push_back(geometry[*vertexIter]);

    return coords;
  }

  static std::map<vertex_type, double> getReferenceQuadrature()
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

  std::map<vertex_type, double> getQuadrature(const mesh_geometry<dimension>& geometry) const
  {
    const std::map<vertex_type, double> referenceWeightings(getReferenceQuadrature());
    const double scaling = getArea(geometry) / 0.5; // 0.5 is area of reference triangle
    std::map<vertex_type, double> weightings;

    for(std::map<vertex_type, double>::const_iterator refIter(referenceWeightings.begin()); refIter!=referenceWeightings.end(); ++refIter)
      weightings[reference_to_physical(geometry, refIter->first)] = refIter->second * scaling;

    return weightings;
  }

  void print(std::ostream& out = std::cout) const
  {
    for(unsigned int i=0; i<vertex_ids.size(); ++i)
      out << "Vertex: " << vertex_ids[i] << std::endl;
  }

  double getArea(const mesh_geometry<dimension>& geometry) const
  {
    const std::vector<vertex_type> vertices(getCoordinates(geometry));
    const double doubleArea = vertices[0][0] * (vertices[1][1] - vertices[2][1]) +
                              vertices[1][0] * (vertices[2][1] - vertices[0][1]) +
                              vertices[2][0] * (vertices[0][1] - vertices[1][1]);
    return doubleArea / 2.0;
  }

  std::set< std::pair<vertex_id, vertex_id> > getFacets() const
  {
    std::set< std::pair<vertex_id, vertex_id> > facets;

    for(unsigned i=0; i<vertex_count; ++i)
      facets.insert(std::make_pair(vertex_ids[i], vertex_ids[(i+1)%vertex_count]));
    
    return facets;
  }

  vertex_type reference_to_physical(const mesh_geometry<dimension>& geometry, const vertex_type& vertex) const
  {
    const std::vector<vertex_type> vertices(getCoordinates(geometry));

    const double xsi = vertex[0];
    const double eta = vertex[1];

    const vertex_type v = (1.0 - xsi - eta) * vertices[0] +
                          xsi * vertices[1] +
                          eta *vertices[2];

    return v;
  }

  bool contains(const mesh_geometry<dimension>& geometry, const vertex_type& v) const
  {
    const std::vector<vertex_type> vertices(getCoordinates(geometry));
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

  std::set< std::set<std::size_t> > getIncidentVertices(MeshTopology& topology, const MeshEntity& cellEntity, std::size_t d) const
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
};

}
#endif
