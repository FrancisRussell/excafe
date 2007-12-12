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

  static vertex_type getVertex(const std::map<vertex_id, vertex_type>& vertexMap, const vertex_id vid)
  {
    const std::map<vertex_id, vertex_type>::const_iterator mapIter = vertexMap.find(vid);
    assert(mapIter != vertexMap.end());
    return mapIter->second;
  }

  std::vector<vertex_type> getCoordinates(const std::map<vertex_id, vertex_type>& vertexMap) const
  {
    const std::vector<vertex_id> vertex_ids(getIndices());
    std::vector< vertex<dimension> > coords;

    for(std::vector<vertex_id>::const_iterator vertexIter(vertex_ids.begin()); vertexIter!=vertex_ids.end(); ++vertexIter)
    {
      coords.push_back(getVertex(vertexMap, *vertexIter));
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

  std::map<vertex_type, double> getQuadrature(const std::map<vertex_id, vertex_type>& vertexMap) const
  {
    const std::map<vertex_type, double> referenceWeightings(getReferenceQuadrature());
    const double scaling = getArea(vertexMap) / 0.5; // 0.5 is area of reference triangle
    std::map<vertex_type, double> weightings;

    for(std::map<vertex_type, double>::const_iterator refIter(referenceWeightings.begin()); refIter!=referenceWeightings.end(); ++refIter)
      weightings[reference_to_physical(vertexMap, refIter->first)] = refIter->second * scaling;

    return weightings;
  }

  void print(std::ostream& out = std::cout) const
  {
    for(unsigned int i=0; i<vertex_ids.size(); ++i)
      out << "Vertex: " << vertex_ids[i] << std::endl;
  }

  double getArea(const std::map<vertex_id, vertex_type>& vertexMap) const
  {
    const std::vector<vertex_type> vertices(getCoordinates(vertexMap));
    const double area = vertices[0][0] * (vertices[1][1] - vertices[2][1]) +
                        vertices[1][0] * (vertices[2][1] - vertices[0][1]) +
                        vertices[2][0] * (vertices[0][1] - vertices[1][1]);
    return area;
  }

  vertex_type reference_to_physical(const std::map<vertex_id, vertex_type>& vertexMap, const vertex_type& vertex) const
  {
    const std::vector<vertex_type> vertices(getCoordinates(vertexMap));

    const double xsi = vertex[0];
    const double eta = vertex[1];

    const double x = (1.0 - xsi - eta) * vertices[0][0] +
                     xsi * vertices[1][0] +
                     eta * vertices[2][0];

    const double y = (1.0 - xsi - eta) * vertices[0][0] +
                     xsi * vertices[1][1] +
                     eta * vertices[2][1];

    return vertex_type(x, y);
  }
};

}
#endif
