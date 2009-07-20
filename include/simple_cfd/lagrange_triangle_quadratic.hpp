#ifndef SIMPLE_CFD_LAGRANGE_TRIANGLE_QUADRATIC_HPP
#define SIMPLE_CFD_LAGRANGE_TRIANGLE_QUADRATIC_HPP

#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <cstddef>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include "utility.hpp"
#include "simple_cfd_fwd.hpp"
#include "finite_element.hpp"
#include "dof_numbering_basic.hpp"
#include "dof_association.hpp"
#include "dof.hpp"

namespace cfd
{

template<unsigned int R>
class LagrangeTriangleQuadratic : public FiniteElement<TriangularCell::dimension>
{
public:
  typedef TriangularCell cell_type;
  typedef Tensor<dimension> value_type;
  typedef Tensor<dimension> gradient_type;
  typedef Tensor<dimension> divergence_type;

  static const std::size_t dimension = cell_type::dimension;
  static const std::size_t rank = R;

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

  std::size_t getRank() const
  {
    return rank;
  }

  std::size_t getDimension() const
  {
    return dimension;
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

  value_type evaluateTensor(const CellVertices<dimension>& cellVertices, const std::size_t i, const vertex_type& vRef) const
  {
    assert(i < spaceDimension());
    const DofAssociation dofAssociation = dofNumbering.getLocalAssociation(i);
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);
    const unsigned node_on_cell = dofAssociation.getEntityDimension() == 0 ? 
      dofAssociation.getEntityIndex() : dofAssociation.getEntityIndex()+3;
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

    value_type result(rank);
    result[tensorIndex.data()] = (gf/gn) * (hf/hn);

    return result;
  }

  gradient_type evaluateGradient(const CellVertices<dimension>& cellVertices, const std::size_t i, const vertex_type& vRef) const
  {
    assert(i < spaceDimension());

    const DofAssociation dofAssociation = dofNumbering.getLocalAssociation(i);
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);
    const unsigned node_on_cell = dofAssociation.getEntityDimension() == 0 ? 
      dofAssociation.getEntityIndex() : dofAssociation.getEntityIndex()+3;
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

    gradient_type result(rank + 1);

    result[xTensorIndex.data()] = ((vertices[j2][1] - vertices[j1][1])/gn) * (hf/hn) +
                (gf/gn) * ((vertices[k2][1] - vertices[k1][1])/hn);

    result[yTensorIndex.data()] = -((vertices[j2][0] - vertices[j1][0])/gn) * (hf/hn) -
                 (gf/gn) * ((vertices[k2][0] - vertices[k1][0])/hn);

    return result;
  }

  divergence_type evaluateDivergence(const CellVertices<dimension>& cellVertices, const std::size_t i, const vertex_type& vRef) const
  {
    assert(i < spaceDimension());
    const DofAssociation dofAssociation = dofNumbering.getLocalAssociation(i);
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);
    const unsigned node_on_cell = dofAssociation.getEntityDimension() == 0 ? 
      dofAssociation.getEntityIndex() : dofAssociation.getEntityIndex()+3;
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

    divergence_type result(rank - 1);

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

  unsigned spaceDimension() const
  {
    return 6 * detail::Power<dimension, rank>::value;
  }

  vertex_type getDofCoordinateGlobal(const Mesh<dimension>& m, const cell_id cid, const std::size_t dof) const
  {
    assert((dof>=0 && dof< 6*detail::Power<dimension, rank>::value));
    const CellVertices<dimension> vertices = m.getCoordinates(cid);
    return referenceCell.referenceToPhysical(vertices, getDofCoordinateLocal(dof));
  }

  vertex_type getDofCoordinateLocal(const std::size_t dof) const
  {
    assert((dof>=0 && dof< 6*detail::Power<dimension, rank>::value));
    const DofAssociation association = dofNumbering.getLocalAssociation(dof);

    if (association.getEntityDimension() == 0)
    {
      return referenceCell.getLocalVertex(association.getEntityIndex());
    }
    else
    {
      const int vid1 = (association.getEntityIndex())%3;
      const int vid2 = (association.getEntityIndex()+1)%3;
      return (referenceCell.getLocalVertex(vid1) + referenceCell.getLocalVertex(vid2))/2.0;
    }
  }

  // NOTE: by permitting mapping dofs to tensor indices, this commits
  // us to using standard bases.
  unsigned getTensorIndex(const Mesh<dimension>& mesh, const std::size_t cid, const std::size_t dof) const
  {
    assert(dof < spaceDimension());
    return dofNumbering.getTensorIndex(dof);
  }

  virtual std::vector< std::set<dof_t> > resolveIdenticalDofs(const Mesh<dimension>& m, const MeshEntity& entity, const std::set<dof_t>& dofsOnEntity) const
  {
    typedef std::map<std::size_t, std::set<dof_t> > tensor_index_to_dofs_map;
    tensor_index_to_dofs_map tensorIndexToDofsMap;

    BOOST_FOREACH(const dof_t& dof, dofsOnEntity)
    {
      assert(dof.getElement() == this);
      tensorIndexToDofsMap[dofNumbering.getTensorIndex(dof.getIndex())].insert(dof);
    }

    std::vector< std::set<dof_t> > sharedDofs;
    BOOST_FOREACH(const tensor_index_to_dofs_map::value_type& indexMapping, tensorIndexToDofsMap)
    {
      sharedDofs.push_back(indexMapping.second);
    }
    return sharedDofs;
  }

  virtual std::set< Dof<dimension> > getDofsOnEntity(MeshTopology& topology, const cell_id cid, const MeshEntity& entity) const
  {
    const std::size_t space_dimension = spaceDimension();
    const std::size_t localIndex = referenceCell.getLocalIndex(topology, cid, entity);
    const MeshEntity localEntity = MeshEntity(entity.getDimension(), localIndex);

    std::set< Dof<dimension> > result;

    for(std::size_t dof=0; dof<space_dimension; ++dof)
    {
      if (dofNumbering.getLocalAssociation(dof).getEntity() == localEntity)
      {
        result.insert(Dof<dimension>(this, cid, dof));
      }
    }

    return result;
  }
};

}

#endif
