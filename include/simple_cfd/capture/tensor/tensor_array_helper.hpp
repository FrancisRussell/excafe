#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_HELPER_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_HELPER_HPP

#include <cstddef>
#include "index.hpp"
#include "index_generator.hpp"
#include "tensor_array_ref.hpp"
#include "tensor_array_collective_sum.hpp"
#include "tensor_array_collective_product.hpp"
#include "tensor_array_table_polynomial.hpp"
#include "tensor_array_placeholder_position.hpp"
#include "tensor_array_placeholder_vertices.hpp"
#include "tensor_array_matrix.hpp"
#include "array_index_summation.hpp"
#include <simple_cfd/exception.hpp>
#include <simple_cfd/cell_manager.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class TensorArrayHelper
{
public:
  static const std::size_t dimension = D;

private:
  typedef typename CellManager::ref<dimension>::general cell_ref_t;
  cell_ref_t cell;
  IndexGenerator generator;
  TensorArrayRef position;
  ArrayIndexVariable cellVertexIndex;
  TensorArrayRef cellVertices;

  TensorArrayRef jacobian;
  TensorArrayRef jacobianDeterminant;

  bool equalSize(const TensorArrayRef& a, const TensorArrayRef& b) const
  {
    return a->getTensorSize() == b->getTensorSize();
  }

  bool equalDimension(const TensorArrayRef& a, const TensorArrayRef& b) const
  {
    return a->getTensorSize().getDimension() == b->getTensorSize().getDimension();
  }

  TensorArrayRef buildJacobian()
  {
    const FiniteElement<dimension>& coordinateMapping = cell->getCoordinateMapping();
    assert(coordinateMapping.getRank() == 0);

    //FIXME: we currently assume space dimension == number of vertices, for a linear mapping.
    const TensorArrayRef interpolationBases = coordinateMapping.getBasisFunctions(generator,
      cellVertexIndex, position);

    const TensorSize vertexSize(1, dimension);
    const TensorSize coeffSize(0, dimension);

    TensorArrayCollectiveProduct scaledVertices(vertexSize);
    const TensorIndex vertexIndex = scaledVertices.getVisibleIndices();
    const TensorIndex coeffIndex(coeffSize);
    scaledVertices.addOperand(coeffIndex, interpolationBases);
    scaledVertices.addOperand(vertexIndex, cellVertices);
    const ArrayIndexSummation interpolatedVertex(TensorArrayRef::cloneFrom(scaledVertices), cellVertexIndex);
    return localGradient(TensorArrayRef::cloneFrom(interpolatedVertex));
  }

  TensorArrayRef asRankZeroTensor(const TensorArray::polynomial_t& s)
  {
    const ArraySize nullArraySize(0);
    const TensorSize tensorSize(0, dimension);
    TensorArrayTablePolynomial scalarTable(generator, nullArraySize, tensorSize);

    const ArrayIndex arrayIndex(nullArraySize);
    const TensorIndex tensorIndex(tensorSize);
    scalarTable(arrayIndex, tensorIndex) = s;

    return TensorArrayRef::cloneFrom(scalarTable);
  }

public:
  TensorArrayHelper(const cell_ref_t _cell) : 
    cell(_cell), position(new TensorArrayPlaceholderPosition(dimension)), 
    cellVertexIndex(generator.newArrayIndexVariable(cell->numEntities(dimension))), 
    cellVertices(new TensorArrayPlaceholderVertices(cellVertexIndex, dimension)),
    jacobian(buildJacobian()), jacobianDeterminant(determinant(jacobian))
  {
  }

  TensorArrayRef add(const TensorArrayRef& a, const TensorArrayRef& b)
  {
    if (!equalSize(a, b)) CFD_EXCEPTION("Tensor addition operands must be same size.");
    const TensorSize resultSize = a->getTensorSize();

    TensorArrayCollectiveSum sum(resultSize);
    const TensorIndex index = sum.getVisibleIndices();
    sum.addOperand(index, a);
    sum.addOperand(index, b);

    return TensorArrayRef::cloneFrom(sum);
  }

  TensorArrayRef inner(const TensorArrayRef& a, const TensorArrayRef& b)
  {
    if (!equalDimension(a, b)) CFD_EXCEPTION("Tensor inner product operands must have same dimension.");

    const std::size_t dimension = a->getTensorSize().getDimension();
    const std::size_t aRank = a->getTensorSize().getRank();
    const std::size_t bRank = b->getTensorSize().getRank();

    const TensorSize resultSize(aRank + bRank - 2, dimension);

    TensorArrayCollectiveProduct product(resultSize);
    const TensorIndexVariable summationIndex = product.newHiddenIndex();

    const TensorIndex aIndex = product.getVisibleIndices().head(aRank-1).append(summationIndex);
    const TensorIndex bIndex = product.getVisibleIndices().tail(bRank-1).prepend(summationIndex);
    product.addOperand(aIndex, a);
    product.addOperand(bIndex, b);

    return TensorArrayRef::cloneFrom(product);
  }

  TensorArrayRef outer(const TensorArrayRef& a, const TensorArrayRef& b)
  {
    if (!equalDimension(a, b)) CFD_EXCEPTION("Tensor outer product operands must have same dimension.");

    const std::size_t dimension = a->getTensorSize().getDimension();
    const std::size_t aRank = a->getTensorSize().getRank();
    const std::size_t bRank = b->getTensorSize().getRank();

    const TensorSize resultSize(aRank + bRank, dimension);

    TensorArrayCollectiveProduct product(resultSize);
    const TensorIndex aIndex = product.getVisibleIndices().head(aRank);
    const TensorIndex bIndex = product.getVisibleIndices().tail(bRank);
    product.addOperand(aIndex, a);
    product.addOperand(bIndex, b);

    return TensorArrayRef::cloneFrom(product);
  }

  TensorArrayRef colon(const TensorArrayRef& a, const TensorArrayRef& b)
  {
    if (!equalSize(a, b)) CFD_EXCEPTION("Tensor colon product operands must have same size.");

    const std::size_t rank = a->getTensorSize().getRank();
    const std::size_t dimension = a->getTensorSize().getDimension();
    const TensorSize resultSize(0, dimension);

    TensorArrayCollectiveProduct product(resultSize);
    product.createNewHiddenIndices(rank);
    const TensorIndex index = product.getHiddenIndices();
    product.addOperand(index, a);
    product.addOperand(index, b);

    return TensorArrayRef::cloneFrom(product);
  }

  TensorArrayRef determinant(const TensorArrayRef& t)
  {
    if (t->getTensorSize().getRank() != 2) CFD_EXCEPTION("Can only take determinant of a rank 2 tensor.");

    TensorArrayMatrixRight matrix(t);
    TensorArray::polynomial_t det(0.0);

    if (matrix.getDimension() == 1)
    {
      det = matrix(0,0);
    }
    else if (matrix.getDimension() == 2)
    {
      det = matrix(0,0) * matrix (1,1) - matrix(0,1) * matrix(1,0);
    }
    else
    {
      CFD_EXCEPTION("Determinant only implemented for 1x1 and 2x2 matrices.");
    }

    return asRankZeroTensor(det);
  }

  TensorArrayRef localGradient(const TensorArrayRef& v)
  {
    if (position->getTensorSize().getRank() !=1) CFD_EXCEPTION("Can only take gradient w.r.t. tensor of rank 1.");

    const std::size_t rank = v->getTensorSize().getRank();
    const std::size_t dimension = v->getTensorSize().getDimension();
    const ArraySize tableArraySize(0);
    const TensorSize tableTensorSize(1, dimension);

    IndexGenerator g;
    TensorArrayTablePolynomial result(g, tableArraySize, tableTensorSize);
    result.appendAdditionalTensorIndices(g, rank);

    for(std::size_t d=0; d<rank; ++d)
    {
      const ArrayIndex arrayIndex(tableArraySize);
      TensorIndex tensorIndex(tableTensorSize);
      tensorIndex[0] = d;

      const ScalarPlaceholder derivative(v->derivative(position[d]), result.getAdditionalIndices());
      result(arrayIndex, tensorIndex) = TensorArray::polynomial_t(derivative); 
    }

    return TensorArrayRef::cloneFrom(result);
  }

  TensorArrayRef globalGradient(const TensorArrayRef& v)
  {
    //FIXME: implement me!
  }
};

}

}
#endif
