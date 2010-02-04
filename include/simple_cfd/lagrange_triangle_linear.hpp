#ifndef SIMPLE_CFD_LAGRANGE_TRIANGLE_LINEAR_HPP
#define SIMPLE_CFD_LAGRANGE_TRIANGLE_LINEAR_HPP

#include <map>
#include <vector>
#include <utility>
#include <cassert>
#include <cmath>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include "simple_cfd_fwd.hpp"
#include "mesh.hpp"
#include "finite_element.hpp"
#include "numeric/tensor.hpp"
#include "dof_numbering_basic.hpp"
#include "capture/array/tensor_array_function_polynomial.hpp"

namespace cfd
{

template<unsigned R>
class LagrangeTriangleLinear : public FiniteElement<TriangularCell::dimension>
{
public:
  typedef TriangularCell cell_type;
  static const std::size_t rank = R;
  static const std::size_t dimension = cell_type::dimension;
  typedef Tensor<dimension> value_type;
  typedef Tensor<dimension> gradient_type;
  typedef Tensor<dimension> divergence_type;
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

  std::size_t getRank() const
  {
    return rank;
  }

  std::size_t getDimension() const
  {
    return dimension;
  }

  value_type evaluateTensor(const CellVertices<dimension>& vertices, const std::size_t i, const vertex_type& vRef) const
  {
    assert(i < spaceDimension());

    const double area = referenceCell.getArea(vertices);

    const vertex_type v = referenceCell.referenceToPhysical(vertices, vRef);
    const DofAssociation dofAssociation = dofNumbering.getLocalAssociation(i);

    assert(dofAssociation.getEntityDimension() == 0);
    const unsigned node_on_cell = dofAssociation.getEntityIndex();
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);

    const int ip1 = (node_on_cell+1) % 3;
    const int ip2 = (node_on_cell+2) % 3;

    boost::array<std::size_t, rank> tensorIndex;
    convert_to_tensor_index(index_into_tensor, tensorIndex.data());

    value_type result(rank);
    result[tensorIndex.data()] = ((vertices[ip2][0] - vertices[ip1][0]) * (v[1] - vertices[ip1][1]) -
                          (vertices[ip2][1] - vertices[ip1][1]) * (v[0] - vertices[ip1][0])) / (2.0 * area);

    return result;
  }

  gradient_type evaluateGradient(const CellVertices<dimension>& vertices, const std::size_t i, const vertex_type& vRef) const
  {
    assert(i < spaceDimension());
    const double area = referenceCell.getArea(vertices);

    const vertex_type v = referenceCell.referenceToPhysical(vertices, vRef);
    const DofAssociation dofAssociation = dofNumbering.getLocalAssociation(i);

    assert(dofAssociation.getEntityDimension() == 0);
    const unsigned node_on_cell = dofAssociation.getEntityIndex();
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);

    const int ip1 = (node_on_cell+1) % 3;
    const int ip2 = (node_on_cell+2) % 3;

    boost::array<std::size_t, rank+1> xTensorIndex, yTensorIndex;
    convert_to_tensor_index(index_into_tensor, xTensorIndex.data()+1);
    convert_to_tensor_index(index_into_tensor, yTensorIndex.data()+1);
    xTensorIndex[0] = 0;
    yTensorIndex[0] = 1;

    gradient_type result(rank + 1);
    result[xTensorIndex.data()] = -(vertices[ip2][1] - vertices[ip1][1]) / (2.0 * area);
    result[yTensorIndex.data()] =  (vertices[ip2][0] - vertices[ip1][0]) / (2.0 * area);

    return result;
  }

  divergence_type evaluateDivergence(const CellVertices<dimension>& vertices, const std::size_t i, const vertex_type& vRef) const
  {
    assert(i < spaceDimension());

    const double area = referenceCell.getArea(vertices);

    const vertex_type v = referenceCell.referenceToPhysical(vertices, vRef);
    const DofAssociation dofAssociation = dofNumbering.getLocalAssociation(i);

    assert(dofAssociation.getEntityDimension() == 0);
    const unsigned node_on_cell = dofAssociation.getEntityIndex();
    const unsigned index_into_tensor = dofNumbering.getTensorIndex(i);

    const int ip1 = (node_on_cell+1) % 3;
    const int ip2 = (node_on_cell+2) % 3;

    // Note how we don't use the final value in the index
    // FIXME: if we don't use an index of the original size, we'll buffer overflow
    boost::array<std::size_t, rank> tensorIndex;
    convert_to_tensor_index(index_into_tensor, tensorIndex.data());

    divergence_type result(rank - 1);

    if (tensorIndex[0] == 0)
      result[tensorIndex.data()+1] += -(vertices[ip2][1] - vertices[ip1][1]) / (2.0 * area);
    else if (tensorIndex[0] == 1)
      result[tensorIndex.data()+1] +=  (vertices[ip2][0] - vertices[ip1][0]) / (2.0 * area);
    else
      assert(false);

    return result;
  }

  unsigned spaceDimension() const
  {
    return 3 * detail::Power<dimension, rank>::value;
  }

  vertex_type getDofCoordinateGlobal(const Mesh<dimension>& m, const cell_id cid, const std::size_t dof) const
  {
    assert(dof>=0 && dof<(3 * detail::Power<dimension, rank>::value));
    const CellVertices<dimension> vertices = m.getCoordinates(cid);
    return referenceCell.referenceToPhysical(vertices, getDofCoordinateLocal(dof));
  }

  vertex_type getDofCoordinateLocal(const std::size_t dof) const
  {
    assert(dof>=0 && dof<(3 * detail::Power<dimension, rank>::value));
    const DofAssociation association = dofNumbering.getLocalAssociation(dof);
    return referenceCell.getLocalVertex(association.getEntityIndex());
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

  virtual const GeneralCell<dimension>& getCell() const
  {
    return referenceCell;
  }
};

}

#endif
