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
#include "dof_numbering_basic.hpp"
#include <cmath>

namespace cfd
{

template<unsigned R>
class LagrangeTriangleLinear : public FiniteElement<TriangularCell>
{
public:
  typedef TriangularCell cell_type;
  static const std::size_t rank = R;
  static const std::size_t dimension = cell_type::dimension;
  typedef Tensor<dimension, rank, double> value_type;
  typedef Tensor<dimension, rank+1, double> gradient_type;
  typedef Tensor<dimension, rank-1, double> divergence_type;
  typedef vertex<dimension> vertex_type;

private:
  static const unsigned int tensor_size = detail::Power<dimension, rank>::value;
  static const unsigned int dofs_per_index = 3;
  const cell_type referenceCell;
  DofNumberingBasic<dimension> dofNumbering;

  // This converts a value to a list of tensor indices in row major order.
  // The order is irrelevent so long as it is consistent and can be used to
  // determine common DoFs between cells
  static void convert_to_tensor_index(const unsigned index, std::size_t* indices)
  {
    unsigned remainder = index;

    for(std::size_t i=0; i<rank; ++i)
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
    dofsPerEntity[1] = 0; 
    dofsPerEntity[2] = 0; 
    return DofNumberingBasic<dimension>(referenceCell, dofsPerEntity, tensor_size);
  }

public:
  // We define the numbering of bases on a cell in the following fashion
  // index_into_tensor * number_of_nodes_on_cell + node_on_cell_id

  LagrangeTriangleLinear() : dofNumbering(buildDofNumberingHelper())
  {
  }

  value_type evaluate_tensor(const CellVertices<dimension>& vertices, const unsigned int i, const vertex_type& vRef) const
  {
    assert(i < space_dimension());

    const double area = referenceCell.getArea(vertices);

    const vertex_type v = referenceCell.referenceToPhysical(vertices, vRef);
    const std::pair<MeshEntity, std::size_t> dofLocation = dofNumbering.getLocalLocation(i);
    const unsigned node_on_cell = dofLocation.first.getIndex();
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);

    const int ip1 = (node_on_cell+1) % 3;
    const int ip2 = (node_on_cell+2) % 3;

    boost::array<std::size_t, rank> tensorIndex;
    convert_to_tensor_index(index_into_tensor, tensorIndex.data());

    value_type result;
    result[tensorIndex.data()] = ((vertices[ip2][0] - vertices[ip1][0]) * (v[1] - vertices[ip1][1]) -
                          (vertices[ip2][1] - vertices[ip1][1]) * (v[0] - vertices[ip1][0])) / (2.0 * area);

    return result;
  }

  gradient_type evaluate_gradient(const CellVertices<dimension>& vertices, const unsigned int i, const vertex_type& vRef) const
  {
    assert(i < space_dimension());
    const double area = referenceCell.getArea(vertices);

    const vertex_type v = referenceCell.referenceToPhysical(vertices, vRef);
    const std::pair<MeshEntity, std::size_t> dofLocation = dofNumbering.getLocalLocation(i);
    const unsigned node_on_cell = dofLocation.first.getIndex();
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);

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

  divergence_type evaluate_divergence(const CellVertices<dimension>& vertices, const unsigned int i, const vertex_type& vRef) const
  {
    assert(i < space_dimension());

    const double area = referenceCell.getArea(vertices);

    const vertex_type v = referenceCell.referenceToPhysical(vertices, vRef);
    const std::pair<MeshEntity, std::size_t> dofLocation = dofNumbering.getLocalLocation(i);
    const unsigned node_on_cell = dofLocation.first.getIndex();
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);

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

  std::vector< std::pair<unsigned, unsigned> > getCommonDegreesOfFreedom(const Mesh<dimension>& m, const cell_id cid, const cell_id cid2) const
  {
    const std::vector<vertex_id> cid_vertices(m.getIndices(MeshEntity(dimension, cid), 0));
    const std::vector<vertex_id> cid2_vertices(m.getIndices(MeshEntity(dimension, cid2), 0));

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

  vertex_type getDofCoordinateGlobal(const Mesh<dimension>& m, const cell_id cid, const unsigned dof) const
  {
    assert(dof>=0 && dof<(3 * detail::Power<dimension, rank>::value));
    const CellVertices<dimension> vertices = m.getCoordinates(cid);
    return referenceCell.referenceToPhysical(vertices, getDofCoordinateLocal(dof));
  }

  vertex_type getDofCoordinateLocal(const unsigned dof) const
  {
    assert(dof>=0 && dof<(3 * detail::Power<dimension, rank>::value));
    const std::pair<MeshEntity, std::size_t> location = dofNumbering.getLocalLocation(dof);
    return referenceCell.getLocalVertex(location.first.getIndex());
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
    const std::size_t localIndex = referenceCell.getLocalIndex(topology, cid, entity);
    std::set< boost::tuple<const FiniteElement<cell_type>*, cell_id, std::size_t> > result;

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
