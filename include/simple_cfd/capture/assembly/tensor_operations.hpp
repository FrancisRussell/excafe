#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_TENSOR_OPERATIONS_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_TENSOR_OPERATIONS_HPP

#include <stack>
#include <algorithm>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/numeric/tensor.hpp>
#include <simple_cfd/numeric/tensor_matrix_view.hpp>
#include <simple_cfd/numeric/functional.hpp>
#include "scalar_placeholder.hpp"
#include "position_placeholder.hpp"
#include "cell_vertices_placeholder.hpp"
#include "scalar_access.hpp"

namespace cfd
{

namespace detail
{

template<std::size_t D>
class TensorOperations
{
private:
  static const std::size_t dimension = D;
  typedef ScalarPlaceholder::expression_t expression_t;
  typedef Tensor<dimension, expression_t> tensor_t;
  PositionPlaceholder position;

public:
  tensor_t grad(const tensor_t& operand) const
  {
    const TensorSize resultSize(operand.getRank()+1, dimension);
    tensor_t result(resultSize);

    for(std::size_t d=0; d<dimension; ++d)
    {
      const ScalarPlaceholder coord(position[d]);
      tensor_t derivative(operand);
      std::transform(derivative.begin(), derivative.end(), derivative.begin(),
        PolynomialDifferentiator<expression_t>(coord));
      result.setElement(d, derivative);
    }

    return result;
  }

  tensor_t adjugate(const tensor_t& t) const
  {
    tensor_t operand(t);
    TensorMatrixView<dimension, expression_t> operandMatrix(operand);

    tensor_t adjugateTensor(operand.getSize());
    TensorMatrixView<dimension, expression_t> adjugateMatrix(adjugateTensor);

    if (operand.getDimension() == 1)
    {
      adjugateMatrix(0,0) = operandMatrix(0,0);
    }
    else if (operand.getDimension() == 2)
    {
      adjugateMatrix(0,0) = operandMatrix(1,1);
      adjugateMatrix(0,1) = operandMatrix(0,1) * -1.0;
      adjugateMatrix(1,0) = operandMatrix(1,0) * -1.0;
      adjugateMatrix(1,1) = operandMatrix(0,0);
    }
    else
    {
      CFD_EXCEPTION("Determinant only implemented for 1x1 and 2x2 matrices.");
    }

    return adjugateTensor; 
  }

  tensor_t transpose(const tensor_t& t) const
  {
    tensor_t operand(t);
    TensorMatrixView<dimension, expression_t> operandMatrix(operand);

    tensor_t transposeTensor(operand.getSize());
    TensorMatrixView<dimension, expression_t> transposeMatrix(transposeTensor);

    for(std::size_t row=0; row<dimension; ++row)
    {
      for(std::size_t col=0; col<dimension; ++col)
      {
        transposeMatrix(col, row) = operandMatrix(row, col);
      }
    }

    return transposeTensor;
  }

  tensor_t invert(const tensor_t& t) const
  {
    return adjugate(t)/determinant(t);
  }

  expression_t determinant(const tensor_t& t) const
  {
    tensor_t operand(t);
    TensorMatrixView<dimension, expression_t> operandMatrix(operand);
    expression_t det(0.0);

    if (operand.getDimension() == 1)
    {
      det = operandMatrix(0,0);
    }
    else if (operand.getDimension() == 2)
    {
      det = operandMatrix(0,0) * operandMatrix (1,1) - operandMatrix(0,1) * operandMatrix(1,0);
    }
    else
    {
      CFD_EXCEPTION("Determinant only implemented for 1x1 and 2x2 matrices.");
    }

    return expression_t::group(det);
  }
  
  tensor_t gradToDiv(const tensor_t& operand) const
  {
    // In order to take div, original operand must have had rank >= 1, so gradient rank must be>= 2.
    assert(operand.getRank() >= 2);
    const TensorSize gradientSize(operand.getSize());
    const TensorSize divergenceSize(gradientSize.getRank() - 2, gradientSize.getDimension());
    tensor_t div(divergenceSize);
    
    for(std::size_t divIndexFlat = 0; divIndexFlat < div.getExtent(); ++divIndexFlat)
    {
      for(std::size_t d=0; d<dimension; ++d)
      {
        // First two indices into gradient are equal
        const std::size_t gradIndexFlat = (d*dimension + d) * divergenceSize.getExtent() + divIndexFlat;
        const TensorIndex divIndex = TensorIndex::unflatten(divergenceSize, divIndexFlat, row_major_tag());
        const TensorIndex gradIndex = TensorIndex::unflatten(gradientSize, gradIndexFlat, row_major_tag());

        div[divIndex] += operand[gradIndex];
      }
    }

    return div;
  }
};

}

}

#endif
