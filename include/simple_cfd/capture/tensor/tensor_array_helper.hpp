#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_HELPER_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_HELPER_HPP

#include "index.hpp"
#include "index_generator.hpp"
#include "tensor_array_ref.hpp"
#include "tensor_array_collective_sum.hpp"
#include "tensor_array_collective_product.hpp"
#include "tensor_array_table_polynomial.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace detail
{

class TensorArrayHelper
{
private:
  bool equalSize(const TensorArrayRef& a, const TensorArrayRef& b) const
  {
    return a->getTensorSize() == b->getTensorSize();
  }

  bool equalDimension(const TensorArrayRef& a, const TensorArrayRef& b) const
  {
    return a->getTensorSize().getDimension() == b->getTensorSize().getDimension();
  }

public:
  TensorArrayRef add(const TensorArrayRef& a, const TensorArrayRef& b) const
  {
    if (!equalSize(a, b)) CFD_EXCEPTION("Tensor addition operands must be same size.");
    const TensorSize resultSize = a->getTensorSize();

    TensorArrayCollectiveSum sum(resultSize);
    const TensorIndex index = sum.getVisibleIndices();
    sum.addOperand(index, a);
    sum.addOperand(index, b);
    return sum;
  }

  TensorArrayRef inner(const TensorArrayRef& a, const TensorArrayRef& b) const
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
    return product;
  }

  TensorArrayRef outer(const TensorArrayRef& a, const TensorArrayRef& b) const
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
    return product;
  }

  TensorArrayRef colon(const TensorArrayRef& a, const TensorArrayRef& b) const
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
    return product;
  }

  TensorArrayRef localGradient(const TensorArrayRef& v, const TensorPlaceholder& pos) const
  {
    if (pos.getTensorSize().getRank() !=1) CFD_EXCEPTION("Can only take gradient w.r.t tensor of rank 1");

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

      const ScalarPlaceholder derivative(v->derivative(pos[d]), result.getAdditionalIndices());
      result(arrayIndex, tensorIndex) = TensorArray::polynomial_t(derivative); 
    }

    return result;
  }
};

}

}
#endif
