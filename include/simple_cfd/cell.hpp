#ifndef SIMPLE_CFD_CELL_HPP
#define SIMPLE_CFD_CELL_HPP

#include<vector>
#include<cassert>
#include<map>
#include<boost/array.hpp>
#include<ostream>
#include<iostream>
#include"mesh.hpp"
#include"vertex.hpp"
#include"simple_cfd_fwd.hpp"

namespace cfd
{

template<shape s>
class cell
{
};

template<>
class cell<triangle>
{
public:
  static const shape cell_shape = triangle;
  static const unsigned int dimension = shape_dimensions<triangle>::dimension;
  typedef vertex<dimension> vertex_type;
  static const unsigned int vertex_count = 3;

private:
  boost::array<vertex_id, vertex_count> vertex_ids;

  std::vector<vertex_type> getCoordinates(const mesh_geometry<dimension>& geometry) const
  {
    const std::vector<vertex_id> vertex_ids(getIndices());
    std::vector< vertex<dimension> > coords;

    for(std::vector<vertex_id>::const_iterator vertexIter(vertex_ids.begin()); vertexIter!=vertex_ids.end(); ++vertexIter)
    {
      coords.push_back(geometry[*vertexIter]);
    }
    return coords;
  }

public:
  cell(const std::vector<vertex_id>& vertices)
  {
    assert(vertices.size() == vertex_count);
    std::copy(vertices.begin(), vertices.end(), vertex_ids.begin());
  }

  std::vector<vertex_id> getIndices() const
  {
    const std::vector<vertex_id> indices(vertex_ids.begin(), vertex_ids.end());
    return indices;
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
};

}
#endif
