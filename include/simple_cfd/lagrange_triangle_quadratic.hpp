#ifndef SIMPLE_CFD_LAGRANGE_TRIANGLE_QUADRATIC_HPP
#define SIMPLE_CFD_LAGRANGE_TRIANGLE_QUADRATIC_HPP

#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <cstddef>
#include <boost/tuple/tuple.hpp>
#include <boost/array.hpp>
#include "utility.hpp"
#include "simple_cfd_fwd.hpp"
#include "finite_element.hpp"
#include "dof_numbering_basic.hpp"

namespace cfd
{

template<unsigned int R>
class LagrangeTriangleQuadratic : public FiniteElement<TriangularCell>
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
  static const unsigned int dofs_per_index = 6;
  const cell_type referenceCell;
  DofNumberingBasic<dimension> dofNumbering;

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

  DofNumberingBasic<dimension> buildDofNumberingHelper() const
  {
    boost::array<std::size_t, dimension+1> dofsPerEntity;
    dofsPerEntity[0] = 1; 
    dofsPerEntity[1] = 1; 
    dofsPerEntity[2] = 0; 
    return DofNumberingBasic<dimension>(referenceCell, dofsPerEntity, tensor_size);
  }

public:
  LagrangeTriangleQuadratic() : dofNumbering(buildDofNumberingHelper())
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

  value_type evaluate_tensor(const CellVertices<dimension>& cellVertices, const unsigned int i, const vertex_type& vRef) const
  {
    assert(i < space_dimension());
    const std::pair<MeshEntity, std::size_t> dofLocation = dofNumbering.getLocalLocation(i);
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);
    const unsigned node_on_cell = dofLocation.first.getDimension() == 0 ? dofLocation.first.getIndex() : dofLocation.first.getIndex()+3;
    const vertex_type v = referenceCell.referenceToPhysical(cellVertices, vRef);

    boost::array<vertex_type, 6> vertices;
    std::copy(cellVertices.begin(), cellVertices.end(), vertices.begin());

    // Create interpolated vertices
    vertices[3] = (vertices[0] + vertices[1])/2;
    vertices[4] = (vertices[1] + vertices[2])/2;
    vertices[5] = (vertices[2] + vertices[0])/2;

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

  gradient_type evaluate_gradient(const CellVertices<dimension>& cellVertices, const unsigned int i, const vertex_type& vRef) const
  {
    assert(i < space_dimension());
    const std::pair<MeshEntity, std::size_t> dofLocation = dofNumbering.getLocalLocation(i);
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);
    const unsigned node_on_cell = dofLocation.first.getDimension() == 0 ? dofLocation.first.getIndex() : dofLocation.first.getIndex()+3;
    const vertex_type v = referenceCell.referenceToPhysical(cellVertices, vRef);

    boost::array<vertex_type, 6> vertices;
    std::copy(cellVertices.begin(), cellVertices.end(), vertices.begin());

    // Create interpolated vertices
    vertices[3] = (vertices[0] + vertices[1])/2;
    vertices[4] = (vertices[1] + vertices[2])/2;
    vertices[5] = (vertices[2] + vertices[0])/2;

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

  divergence_type evaluate_divergence(const CellVertices<dimension>& cellVertices, const unsigned int i, const vertex_type& vRef) const
  {
    assert(i < space_dimension());
    const std::pair<MeshEntity, std::size_t> dofLocation = dofNumbering.getLocalLocation(i);
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);
    const unsigned node_on_cell = dofLocation.first.getDimension() == 0 ? dofLocation.first.getIndex() : dofLocation.first.getIndex()+3;
    const vertex_type v = referenceCell.referenceToPhysical(cellVertices, vRef);

    boost::array<vertex_type, 6> vertices;
    std::copy(cellVertices.begin(), cellVertices.end(), vertices.begin());

    // Create interpolated vertices
    vertices[3] = (vertices[0] + vertices[1])/2;
    vertices[4] = (vertices[1] + vertices[2])/2;
    vertices[5] = (vertices[2] + vertices[0])/2;

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

  std::vector< std::pair<unsigned, unsigned> > getCommonDegreesOfFreedom(const Mesh<dimension>& m, const cell_id cid, const cell_id cid2) const
  {
    const std::vector<vertex_id> cid_vertices(m.getIndices(MeshEntity(dimension, cid), 0));
    const std::vector<vertex_id> cid2_vertices(m.getIndices(MeshEntity(dimension, cid2), 0));

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

  vertex_type getDofCoordinateGlobal(const Mesh<dimension>& m, const cell_id cid, const unsigned dof) const
  {
    assert((dof>=0 && dof< 6*detail::Power<dimension, rank>::value));
    const CellVertices<dimension> vertices = m.getCoordinates(cid);
    return referenceCell.referenceToPhysical(vertices, getDofCoordinateLocal(dof));
  }

  vertex_type getDofCoordinateLocal(const unsigned dof) const
  {
    assert((dof>=0 && dof< 6*detail::Power<dimension, rank>::value));

    if (dof%6 < 3)
    {
      return referenceCell.getLocalVertex(dof%6);
    }
    else
    {
      const int vid1 = (dof%6 - 3)%3;
      const int vid2 = (dof%6 - 2)%3;
      return (referenceCell.getLocalVertex(vid1) + referenceCell.getLocalVertex(vid2))/2.0;
    }
  }

  // NOTE: by permitting mapping dofs to tensor indices, this commits
  // us to using standard bases.
  unsigned getTensorIndex(const unsigned dof) const
  {
    assert(dof < space_dimension());
    return dofNumbering.getTensorIndex(dof);
  }

  virtual std::set< boost::tuple<const FiniteElement<cell_type>*, cell_id, std::size_t> > getDegreesOfFreedom(MeshTopology& topology, const cell_id cid, const MeshEntity& entity) const
  {
    const std::size_t entityIndex = referenceCell.getLocalIndex(topology, cid, entity);
    std::set< boost::tuple<const FiniteElement<cell_type>*, cell_id, std::size_t> > result;

    if (entity.getDimension() == 2) return result;

    // If we are looking at edges, indices are incremented by 3
    const std::size_t localIndex = (entity.getDimension() == 0) ? entityIndex : entityIndex+3;
    for(std::size_t index=0; index < tensor_size; ++index)
    {
      result.insert(boost::make_tuple(this, cid, dofs_per_index*index + localIndex));
    }

    return result;
  }

};

}

#endif
