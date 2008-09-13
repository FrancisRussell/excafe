#ifndef SIMPLE_CFD_LAGRANGE_TRIANGLE_QUADRATIC_HPP
#define SIMPLE_CFD_LAGRANGE_TRIANGLE_QUADRATIC_HPP

#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include "utility.hpp"
#include "simple_cfd_fwd.hpp"
#include "finite_element.hpp"

namespace cfd
{

template<unsigned int R>
class lagrange_triangle_quadratic : public finite_element< cell<triangle> >
{
public:
  typedef cell<triangle> cell_type;
  static const unsigned int rank = R;
  static const unsigned int dimension = cell_type::dimension;
  typedef Tensor<dimension, rank, double> value_type;
  typedef Tensor<dimension, rank+1, double> gradient_type;
  typedef Tensor<dimension, rank-1, double> divergence_type;
  typedef vertex<dimension> vertex_type;

private:
  const mesh<cell_type>* m;

  // This converts a value to a list of tensor indices in row major order.
  // The order is irrelevent so long as it is consistent and can be used to
  // determine common DoFs between cells
  static void convert_to_tensor_index(const unsigned index, std::size_t* indices)
  {
    unsigned remainder = index;

    for(unsigned i=0; i<rank; ++i)
    {
      indices[rank-i-1] = remainder % dimension;
      remainder /= dimension;
    }

    // A fail here means the index was too large
    assert(remainder == 0);
  }

public:
  lagrange_triangle_quadratic(const mesh<cell_type>& _m) : m(&_m)
  {
  }

  /*   Triangle Node Numbering for quadratic basis

       LL        

       2

       |\ 
       | \ 
     5 |  \ 4
       |   \ 
       |    \ 
       ------ 
       0  3  1
  */

  evaluated_basis evaluate_basis(const cell_type& c, const unsigned int i, const vertex_type& v) const
  {
    std::vector<vertex_type> vertices(c.getCoordinates(m->getGeometry()));

    // Create interpolated vertices
    vertices.push_back((vertices[0] + vertices[1])/2);
    vertices.push_back((vertices[1] + vertices[2])/2);
    vertices.push_back((vertices[2] + vertices[0])/2);

    int j1, j2, k1, k2;
    if (i < 3)
    {
      j1 = (i+1)%3;
      j2 = (i+2)%3;
      k1 = 3 + i;
      k2 = 3 + (i + 5)%3;
    }
    else
    {
      j1 = i-3;
      j2 = (i-3+2)%3;
      k1 = (i-3+1)%3;
      k2 = (i-3+2)%3;
    }

    const double gf = (v[0] - vertices[j1][0]) * (vertices[j2][1] - vertices[j1][1]) -
                      (vertices[j2][0] - vertices[j1][0]) * (v[1] - vertices[j1][1]);

    const double gn = (vertices[i][0] - vertices[j1][0]) * (vertices[j2][1] - vertices[j1][1]) -
                      (vertices[j2][0] - vertices[j1][0]) * (vertices[i][1] - vertices[j1][1]);

    const double hf = (v[0] - vertices[k1][0]) * (vertices[k2][1] - vertices[k1][1]) -
                      (vertices[k2][0] - vertices[k1][0]) * (v[1] - vertices[k1][1]);

    const double hn = (vertices[i][0] - vertices[k1][0]) * (vertices[k2][1] - vertices[k1][1]) -
                      (vertices[k2][0] - vertices[k1][0]) * (vertices[i][1] - vertices[k1][1]);

    evaluated_basis result;

    result.value = (gf/gn) * (hf/hn);

    result.dx = ((vertices[j2][1] - vertices[j1][1])/gn) * (hf/hn) +
                (gf/gn) * ((vertices[k2][1] - vertices[k1][1])/hn);

    result.dy = -((vertices[j2][0] - vertices[j1][0])/gn) * (hf/hn) -
                 (gf/gn) * ((vertices[k2][0] - vertices[k1][0])/hn);

    return result;
  }

  value_type evaluate_tensor(const cell_type& c, const unsigned int i, const vertex_type& v) const
  {
    assert(i < space_dimension());
    const unsigned node_on_cell = i % 6;
    const unsigned index_into_tensor = i / 6;
    std::vector<vertex_type> vertices(c.getCoordinates(m->getGeometry()));

    // Create interpolated vertices
    vertices.push_back((vertices[0] + vertices[1])/2);
    vertices.push_back((vertices[1] + vertices[2])/2);
    vertices.push_back((vertices[2] + vertices[0])/2);

    int j1, j2, k1, k2;
    if (node_on_cell < 3)
    {
      j1 = (node_on_cell+1)%3;
      j2 = (node_on_cell+2)%3;
      k1 = 3 + node_on_cell;
      k2 = 3 + (node_on_cell + 5)%3;
    }
    else
    {
      j1 = node_on_cell-3;
      j2 = (node_on_cell-3+2)%3;
      k1 = (node_on_cell-3+1)%3;
      k2 = (node_on_cell-3+2)%3;
    }

    const double gf = (v[0] - vertices[j1][0]) * (vertices[j2][1] - vertices[j1][1]) -
                      (vertices[j2][0] - vertices[j1][0]) * (v[1] - vertices[j1][1]);

    const double gn = (vertices[node_on_cell][0] - vertices[j1][0]) * (vertices[j2][1] - vertices[j1][1]) -
                      (vertices[j2][0] - vertices[j1][0]) * (vertices[node_on_cell][1] - vertices[j1][1]);

    const double hf = (v[0] - vertices[k1][0]) * (vertices[k2][1] - vertices[k1][1]) -
                      (vertices[k2][0] - vertices[k1][0]) * (v[1] - vertices[k1][1]);

    const double hn = (vertices[node_on_cell][0] - vertices[k1][0]) * (vertices[k2][1] - vertices[k1][1]) -
                      (vertices[k2][0] - vertices[k1][0]) * (vertices[node_on_cell][1] - vertices[k1][1]);

    boost::array<std::size_t, rank> tensorIndex;
    convert_to_tensor_index(index_into_tensor, tensorIndex.data());

    value_type result;
    result[tensorIndex.data()] = (gf/gn) * (hf/hn);

    return result;
  }

  gradient_type evaluate_gradient(const cell_type& c, const unsigned int i, const vertex_type& v) const
  {
    assert(i < space_dimension());
    const unsigned node_on_cell = i % 6;
    const unsigned index_into_tensor = i / 6;
    std::vector<vertex_type> vertices(c.getCoordinates(m->getGeometry()));

    // Create interpolated vertices
    vertices.push_back((vertices[0] + vertices[1])/2);
    vertices.push_back((vertices[1] + vertices[2])/2);
    vertices.push_back((vertices[2] + vertices[0])/2);

    int j1, j2, k1, k2;
    if (node_on_cell < 3)
    {
      j1 = (node_on_cell+1)%3;
      j2 = (node_on_cell+2)%3;
      k1 = 3 + node_on_cell;
      k2 = 3 + (node_on_cell + 5)%3;
    }
    else
    {
      j1 = node_on_cell-3;
      j2 = (node_on_cell-3+2)%3;
      k1 = (node_on_cell-3+1)%3;
      k2 = (node_on_cell-3+2)%3;
    }

    const double gf = (v[0] - vertices[j1][0]) * (vertices[j2][1] - vertices[j1][1]) -
                      (vertices[j2][0] - vertices[j1][0]) * (v[1] - vertices[j1][1]);

    const double gn = (vertices[node_on_cell][0] - vertices[j1][0]) * (vertices[j2][1] - vertices[j1][1]) -
                      (vertices[j2][0] - vertices[j1][0]) * (vertices[node_on_cell][1] - vertices[j1][1]);

    const double hf = (v[0] - vertices[k1][0]) * (vertices[k2][1] - vertices[k1][1]) -
                      (vertices[k2][0] - vertices[k1][0]) * (v[1] - vertices[k1][1]);

    const double hn = (vertices[node_on_cell][0] - vertices[k1][0]) * (vertices[k2][1] - vertices[k1][1]) -
                      (vertices[k2][0] - vertices[k1][0]) * (vertices[node_on_cell][1] - vertices[k1][1]);

    boost::array<std::size_t, rank+1> xTensorIndex, yTensorIndex;
    convert_to_tensor_index(index_into_tensor, xTensorIndex.data()+1);
    convert_to_tensor_index(index_into_tensor, yTensorIndex.data()+1);
    xTensorIndex[0] = 0;
    yTensorIndex[0] = 1;

    gradient_type result;
    result[xTensorIndex.data()] = ((vertices[j2][1] - vertices[j1][1])/gn) * (hf/hn) +
                (gf/gn) * ((vertices[k2][1] - vertices[k1][1])/hn);

    result[yTensorIndex.data()] = -((vertices[j2][0] - vertices[j1][0])/gn) * (hf/hn) -
                 (gf/gn) * ((vertices[k2][0] - vertices[k1][0])/hn);

    return result;
  }

  divergence_type evaluate_divergence(const cell_type& c, const unsigned int i, const vertex_type& v) const
  {
    assert(i < space_dimension());
    const unsigned node_on_cell = i % 6;
    const unsigned index_into_tensor = i / 6;

    std::vector<vertex_type> vertices(c.getCoordinates(m->getGeometry()));

    // Create interpolated vertices
    vertices.push_back((vertices[0] + vertices[1])/2);
    vertices.push_back((vertices[1] + vertices[2])/2);
    vertices.push_back((vertices[2] + vertices[0])/2);

    int j1, j2, k1, k2;
    if (node_on_cell < 3)
    {
      j1 = (node_on_cell+1)%3;
      j2 = (node_on_cell+2)%3;
      k1 = 3 + node_on_cell;
      k2 = 3 + (node_on_cell + 5)%3;
    }
    else
    {
      j1 = node_on_cell-3;
      j2 = (node_on_cell-3+2)%3;
      k1 = (node_on_cell-3+1)%3;
      k2 = (node_on_cell-3+2)%3;
    }

    const double gf = (v[0] - vertices[j1][0]) * (vertices[j2][1] - vertices[j1][1]) -
                      (vertices[j2][0] - vertices[j1][0]) * (v[1] - vertices[j1][1]);

    const double gn = (vertices[node_on_cell][0] - vertices[j1][0]) * (vertices[j2][1] - vertices[j1][1]) -
                      (vertices[j2][0] - vertices[j1][0]) * (vertices[node_on_cell][1] - vertices[j1][1]);

    const double hf = (v[0] - vertices[k1][0]) * (vertices[k2][1] - vertices[k1][1]) -
                      (vertices[k2][0] - vertices[k1][0]) * (v[1] - vertices[k1][1]);

    const double hn = (vertices[node_on_cell][0] - vertices[k1][0]) * (vertices[k2][1] - vertices[k1][1]) -
                      (vertices[k2][0] - vertices[k1][0]) * (vertices[node_on_cell][1] - vertices[k1][1]);

    // Note how we don't use the final value in the index
    // FIXME: if we don't use an index of the original size, we'll buffer overflow
    boost::array<std::size_t, rank> tensorIndex;
    convert_to_tensor_index(index_into_tensor, tensorIndex.data());

    divergence_type result;

    if (tensorIndex[0] == 0)
      result[tensorIndex.data()+1] += ((vertices[j2][1] - vertices[j1][1])/gn) * (hf/hn) +
                                      (gf/gn) * ((vertices[k2][1] - vertices[k1][1])/hn);
    else if (tensorIndex[0] == 1)
      result[tensorIndex.data()+1] += -((vertices[j2][0] - vertices[j1][0])/gn) * (hf/hn) -
                                       (gf/gn) * ((vertices[k2][0] - vertices[k1][0])/hn);
    else
      assert(false);

    return result;
  }

  unsigned space_dimension() const
  {
    return 6 * detail::Power<dimension, rank>::value;
  }

  std::vector< std::pair<unsigned, unsigned> > getCommonDegreesOfFreedom(const cell_id cid, const cell_id cid2) const
  {
    const std::vector<vertex_id> cid_vertices(m->getCell(cid).getIndices());
    const std::vector<vertex_id> cid2_vertices(m->getCell(cid2).getIndices());

    std::vector< std::pair<vertex_id, vertex_id> > cid_edges;
    std::vector< std::pair<vertex_id, vertex_id> > cid2_edges;

    for(unsigned edge=0; edge<3; ++edge)
    {
      cid_edges.push_back(std::make_pair(cid_vertices[edge], cid_vertices[(edge+1)%3]));
      cid2_edges.push_back(std::make_pair(cid2_vertices[edge], cid2_vertices[(edge+1)%3]));
    }

    // Map vertices and midpoints of cid2 onto degrees of freedom
    unsigned cid2_dof=0;

    std::map<vertex_id, unsigned> cid2_vertex_dof;
    for(unsigned v=0; v<cid2_vertices.size(); ++v)
    {
      cid2_vertex_dof[cid2_vertices[v]] = cid2_dof;
      ++cid2_dof;
    }

    std::map<std::pair<vertex_id, vertex_id>, unsigned, unordered_pair_compare<vertex_id> > cid2_midpoint_dof;
    for(unsigned m=0; m<cid_edges.size(); ++m)
    {
      cid2_midpoint_dof[cid2_edges[m]] = cid2_dof;
      ++cid2_dof;
    }

    // Iterate over degrees of freedom on cid and find ones that correspond to common vertices and midpoints
    std::vector< std::pair<unsigned, unsigned> > common;
    unsigned cid_dof = 0;

    for(unsigned v=0; v<cid_vertices.size(); ++v)
    {
      const std::map<vertex_id, unsigned>::const_iterator sharedVertexIter = cid2_vertex_dof.find(cid_vertices[v]);

      if (sharedVertexIter != cid2_vertex_dof.end())
        for(unsigned int index_into_tensor = 0; index_into_tensor < detail::Power<dimension, rank>::value; ++index_into_tensor) 
          common.push_back(std::make_pair(index_into_tensor*6 + cid_dof, index_into_tensor*6 + sharedVertexIter->second));

      ++cid_dof;
    }

    for(unsigned m=0; m<cid_edges.size(); ++m)
    {
      const std::map<std::pair<vertex_id, vertex_id>, unsigned>::const_iterator sharedEdgeIter = cid2_midpoint_dof.find(cid_edges[m]);

      if (sharedEdgeIter != cid2_midpoint_dof.end())
        for(unsigned int index_into_tensor = 0; index_into_tensor < detail::Power<dimension, rank>::value; ++index_into_tensor) 
          common.push_back(std::make_pair(index_into_tensor*6 + cid_dof, index_into_tensor*6 + sharedEdgeIter->second));

      ++cid_dof;
    }
    return common;
  }


  std::vector<unsigned> getBoundaryDegreesOfFreedom(const cell_id cid, const std::vector< std::pair<vertex_id, vertex_id> >& boundary) const
  {
    // Create sets for all edges and vertices on boundary
    std::set< std::pair<vertex_id, vertex_id>, unordered_pair_compare<vertex_id> > boundaryEdgeSet(boundary.begin(), boundary.end());
    std::set<vertex_id> boundaryVertices;
    for(std::vector< std::pair<vertex_id, vertex_id> >::const_iterator edgeIter(boundary.begin()); edgeIter!=boundary.end(); ++edgeIter)
    {
      boundaryVertices.insert(edgeIter->first);
      boundaryVertices.insert(edgeIter->second);
    }

    // Create lists of all local edges and vertices
    const std::vector<vertex_id> vertexIndices(m->getCell(cid).getIndices());
    std::vector< std::pair<unsigned, unsigned> > edges;
    for(unsigned edge=0; edge<3; ++edge)
    {
      edges.push_back(std::make_pair(vertexIndices[edge], vertexIndices[(edge+1)%3]));
    }

    // Find degrees of freedom on boundary
    std::vector<unsigned> dofs;
    for(unsigned i=0; i<vertexIndices.size(); ++i)
    {
      if (boundaryVertices.find(vertexIndices[i]) != boundaryVertices.end())
        for(unsigned int index_into_tensor = 0; index_into_tensor < detail::Power<dimension, rank>::value; ++index_into_tensor) 
          dofs.push_back(index_into_tensor*6 + i);
    }

    for(unsigned i=0; i<edges.size(); ++i)
    {
      if (boundaryEdgeSet.find(edges[i]) != boundaryEdgeSet.end())
        for(unsigned int index_into_tensor = 0; index_into_tensor < detail::Power<dimension, rank>::value; ++index_into_tensor) 
          dofs.push_back(index_into_tensor*6 + i+3);
    }

    return dofs;
  }
  
  vertex_type getDofCoordinate(const cell_id cid, const unsigned dof) const
  {
    assert((dof>=0 && dof< 6*detail::Power<dimension, rank>::value));
    if (dof%6 < 3)
    {
      return m->getCoordinates(cid)[dof%6];
    }
    else
    {
      const std::vector<vertex_type> coordinates(m->getCoordinates(cid));
      const int vid1 = (dof%6 - 3)%3;
      const int vid2 = (dof%6 - 2)%3;
      return (coordinates[vid1] + coordinates[vid2])/2;
    }
  }

  // NOTE: by permitting mapping dofs to tensor indices, this commits
  // us to using standard bases.
  unsigned getTensorIndex(const cell_id cid, const unsigned dof) const
  {
    assert(dof < space_dimension());
    return dof/6;
  }
};

}

#endif
