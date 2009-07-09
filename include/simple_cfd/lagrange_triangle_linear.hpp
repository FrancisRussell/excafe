#ifndef SIMPLE_CFD_LAGRANGE_TRIANGLE_LINEAR_HPP
#define SIMPLE_CFD_LAGRANGE_TRIANGLE_LINEAR_HPP

#include <map>
#include <vector>
#include <utility>
#include <cassert>
#include <boost/array.hpp>
#include "simple_cfd_fwd.hpp"
#include "mesh.hpp"
#include "finite_element.hpp"
#include "numeric/tensor.hpp"
#include <cmath>

namespace cfd
{

template<unsigned R>
class lagrange_triangle_linear : public finite_element<TriangularCell>
{
public:
  typedef TriangularCell cell_type;
  static const unsigned int rank = R;
  static const unsigned int dimension = cell_type::dimension;
  typedef Tensor<dimension, rank, double> value_type;
  typedef Tensor<dimension, rank+1, double> gradient_type;
  typedef Tensor<dimension, rank-1, double> divergence_type;
  typedef vertex<dimension> vertex_type;

private:
  static const unsigned int tensor_size = detail::Power<dimension, rank>::value;
  static const unsigned int dofs_per_index = 3;
  const mesh<cell_type>* m;
  const cell_type referenceCell;

  // This converts a value to a list of tensor indices in row major order.
  // The order is irrelevent so long as it is consistent and can be used to
  // determine common DoFs between cells
  static void convert_to_tensor_index(const unsigned index, std::size_t* indices)
  {
    unsigned remainder = index;

    for(int i=0; i<rank; ++i)
    {
      indices[rank-i-1] = remainder % dimension;
      remainder /= dimension;
    }

    // A fail here means the index was too large
    assert(remainder == 0);
  }

public:
  // We define the numbering of bases on a cell in the following fashion
  // index_into_tensor * number_of_nodes_on_cell + node_on_cell_id

  lagrange_triangle_linear(const mesh<cell_type>& _m) : m(&_m)
  {
  }

  value_type evaluate_tensor(const std::size_t cid, const unsigned int i, const vertex_type& vRef) const
  {
    assert(i < space_dimension());
    const vertex_type v = referenceCell.reference_to_physical(*m, cid, vRef);
    const unsigned node_on_cell = i % 3;
    const unsigned index_into_tensor = i / 3;

    const std::vector<vertex_type> vertices(m->getCoordinates(cid));
    const double area = m->getArea(cid);

    const int ip1 = (node_on_cell+1) % 3;
    const int ip2 = (node_on_cell+2) % 3;

    boost::array<std::size_t, rank> tensorIndex;
    convert_to_tensor_index(index_into_tensor, tensorIndex.data());

    value_type result;
    result[tensorIndex.data()] = ((vertices[ip2][0] - vertices[ip1][0]) * (v[1] - vertices[ip1][1]) -
                          (vertices[ip2][1] - vertices[ip1][1]) * (v[0] - vertices[ip1][0])) / (2.0 * area);

    return result;
  }

  gradient_type evaluate_gradient(const std::size_t cid, const unsigned int i, const vertex_type& vRef) const
  {
    assert(i < space_dimension());
    const vertex_type v = referenceCell.reference_to_physical(*m, cid, vRef);
    const unsigned node_on_cell = i % 3;
    const unsigned index_into_tensor = i / 3;

    const std::vector<vertex_type> vertices(m->getCoordinates(cid));
    const double area = m->getArea(cid);

    const int ip1 = (node_on_cell+1) % 3;
    const int ip2 = (node_on_cell+2) % 3;

    boost::array<std::size_t, rank+1> xTensorIndex, yTensorIndex;
    convert_to_tensor_index(index_into_tensor, xTensorIndex.data()+1);
    convert_to_tensor_index(index_into_tensor, yTensorIndex.data()+1);
    xTensorIndex[0] = 0;
    yTensorIndex[0] = 1;

    gradient_type result;
    result[xTensorIndex.data()] = -(vertices[ip2][1] - vertices[ip1][1]) / (2.0 * area);
    result[yTensorIndex.data()] =  (vertices[ip2][0] - vertices[ip1][0]) / (2.0 * area);

    return result;
  }

  divergence_type evaluate_divergence(const std::size_t cid, const unsigned int i, const vertex_type& vRef) const
  {
    assert(i < space_dimension());
    const vertex_type v = referenceCell.reference_to_physical(*m, cid, vRef);
    const unsigned node_on_cell = i % 3;
    const unsigned index_into_tensor = i / 3;

    const std::vector<vertex_type> vertices(m->getCoordinates(cid));
    const double area = m->getArea(cid);

    const int ip1 = (node_on_cell+1) % 3;
    const int ip2 = (node_on_cell+2) % 3;

    // Note how we don't use the final value in the index
    // FIXME: if we don't use an index of the original size, we'll buffer overflow
    boost::array<std::size_t, rank> tensorIndex;
    convert_to_tensor_index(index_into_tensor, tensorIndex.data());

    divergence_type result;

    if (tensorIndex[0] == 0)
      result[tensorIndex.data()+1] += -(vertices[ip2][1] - vertices[ip1][1]) / (2.0 * area);
    else if (tensorIndex[0] == 1)
      result[tensorIndex.data()+1] +=  (vertices[ip2][0] - vertices[ip1][0]) / (2.0 * area);
    else
      assert(false);

    return result;
  }

  unsigned space_dimension() const
  {
    return 3 * detail::Power<dimension, rank>::value;
  }

  std::vector< std::pair<unsigned, unsigned> > getCommonDegreesOfFreedom(const cell_id cid, const cell_id cid2) const
  {
    const std::vector<vertex_id> cid_vertices(m->getIndices(MeshEntity(dimension, cid), 0));
    const std::vector<vertex_id> cid2_vertices(m->getIndices(MeshEntity(dimension, cid2), 0));

    // Map vertices of cid2 onto degrees of freedom
    std::map<vertex_id, unsigned> cid2_dof;
    for(unsigned dof=0; dof<cid2_vertices.size(); ++dof)
      cid2_dof[cid2_vertices[dof]] = dof;

    std::vector< std::pair<unsigned, unsigned> > common;
    // Iterate over degrees of freedom on cid and find ones that correspond to common vertices
    for(unsigned dof=0; dof<cid_vertices.size(); ++dof)
    {
      const std::map<vertex_id, unsigned>::const_iterator sharedVertexIter = cid2_dof.find(cid_vertices[dof]);

      if (sharedVertexIter != cid2_dof.end())
      {
        for(unsigned int index_into_tensor = 0; index_into_tensor < detail::Power<dimension, rank>::value; ++index_into_tensor) 
          common.push_back(std::make_pair(index_into_tensor*3 + dof, index_into_tensor*3 + sharedVertexIter->second));
      }
    }
    return common;
  }

  std::vector<unsigned> getBoundaryDegreesOfFreedom(const cell_id cid, const std::vector< std::pair<vertex_id, vertex_id> >& boundary) const
  {
    // Create set of vertices on boundary
    std::set<vertex_id> boundaryVertices;
    for(std::vector< std::pair<vertex_id, vertex_id> >::const_iterator edgeIter(boundary.begin()); edgeIter!=boundary.end(); ++edgeIter)
    {
      boundaryVertices.insert(edgeIter->first);
      boundaryVertices.insert(edgeIter->second);
    }

    // Create list of local vertices
    const std::vector<vertex_id> vertexIndices(m->getIndices(MeshEntity(dimension, cid), 0));

    // Find vertices on boundary
    std::vector<unsigned> dofs;
    for(unsigned i=0; i<vertexIndices.size(); ++i)
    {
      if (boundaryVertices.find(vertexIndices[i]) != boundaryVertices.end())
      {
        for(unsigned int index_into_tensor = 0; index_into_tensor < detail::Power<dimension, rank>::value; ++index_into_tensor) 
          dofs.push_back(index_into_tensor*3 + i);
      }
    }
    return dofs;
  }

  vertex_type getDofCoordinateGlobal(const cell_id cid, const unsigned dof) const
  {
    assert(dof>=0 && dof<(3 * detail::Power<dimension, rank>::value));
    return referenceCell.reference_to_physical(*m, cid, getDofCoordinateLocal(dof));
  }

  vertex_type getDofCoordinateLocal(const unsigned dof) const
  {
    assert(dof>=0 && dof<(3 * detail::Power<dimension, rank>::value));
    return referenceCell.getLocalVertex(dof % 3);
  }

  // NOTE: by permitting mapping dofs to tensor indices, this commits
  // us to using standard bases.
  unsigned getTensorIndex(const cell_id cid, const unsigned dof) const
  {
    assert(dof < space_dimension());
    return dof/3;
  }

  virtual std::set< boost::tuple<const finite_element<cell_type>*, cell_id, std::size_t> > getDegreesOfFreedom(MeshTopology& topology, const cell_id cid, const MeshEntity& entity) const
  {
    const std::size_t localIndex = referenceCell.getLocalIndex(topology, cid, entity);
    std::set< boost::tuple<const finite_element<cell_type>*, cell_id, std::size_t> > result;

    if (entity.getDimension() == 2 || entity.getDimension() == 1) return result;

    if (entity.getDimension() == 0)
    {
      for(std::size_t index=0; index < tensor_size; ++index)
      {
        result.insert(boost::make_tuple(this, cid, dofs_per_index*index + localIndex));
      }
    }

    return result;
  }
};

}

#endif
