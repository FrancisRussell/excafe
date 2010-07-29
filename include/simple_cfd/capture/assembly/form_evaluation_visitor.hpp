#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_FORM_EVALUATION_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_FORM_EVALUATION_VISITOR_HPP

#include <stack>
#include <algorithm>
#include <simple_cfd/cell_manager.hpp>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/numeric/tensor.hpp>
#include <simple_cfd/numeric/tensor_matrix_view.hpp>
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/numeric/functional.hpp>
#include <simple_cfd/capture/forms/linear_form.hpp>
#include <simple_cfd/capture/forms/field_visitor.hpp>
#include <simple_cfd/capture/forms/field_addition.hpp>
#include <simple_cfd/capture/forms/field_inner_product.hpp>
#include <simple_cfd/capture/forms/field_outer_product.hpp>
#include <simple_cfd/capture/forms/field_colon_product.hpp>
#include <simple_cfd/capture/forms/field_gradient.hpp>
#include <simple_cfd/capture/forms/field_divergence.hpp>
#include <simple_cfd/capture/forms/field_basis.hpp>
#include <simple_cfd/capture/forms/field_discrete_reference.hpp>
#include <simple_cfd/capture/forms/field_scalar.hpp>
#include <simple_cfd/capture/forms/facet_normal.hpp>
#include "scalar_placeholder.hpp"
#include "position_placeholder.hpp"
#include "cell_vertices_placeholder.hpp"
#include "scalar_access.hpp"

namespace cfd
{

namespace detail
{

template<std::size_t D>
class FormEvaluationVisitor : public detail::FieldVisitor
{
private:
  static const std::size_t dimension = D;
  typedef typename CellManager::ref<dimension>::general cell_ref_t;
  typedef ScalarPlaceholder::expression_t expression_t;
  typedef Tensor<dimension, expression_t> tensor_t;

  const Scenario<dimension>& scenario;
  PositionPlaceholder position;
  CellVerticesPlaceholder<dimension> cellVertices;
  const std::size_t basisFunctionIndex;
  std::stack<tensor_t> valueStack;

  tensor_t buildLocalGradient(const tensor_t& operand) const
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

  tensor_t buildGlobalGradient(const tensor_t& operand) const
  {
    const tensor_t localGradient(buildLocalGradient(operand));
    const tensor_t inverseGradient(invert(buildLocalGradient(buildGlobalPosition())));

    return inverseGradient.inner_product(localGradient);
  }

  tensor_t buildGlobalPosition() const
  {
    const cell_ref_t cell = scenario.getMesh().getReferenceCell();
    const std::size_t numVertices = cell->numEntities(0);
    assert(cell->getCoordinateMapping().spaceDimension() == numVertices);

    const TensorSize positionSize(1, dimension);
    tensor_t globalPosition(positionSize);

    for(std::size_t i=0; i<numVertices; ++i)
    {
      globalPosition += cell->getCoordinateMapping().getBasis(i, position) * cellVertices[i];
    }

    return globalPosition;
  }

  tensor_t buildJacobian() const
  {
    return transpose(buildLocalGradient(buildGlobalPosition()));
  }

  static tensor_t adjugate(const tensor_t& t)
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

  static tensor_t transpose(const tensor_t& t)
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

  static tensor_t asTensor(const expression_t& p)
  {
    const TensorSize scalarSize(0, dimension);
    const TensorIndex nullIndex(scalarSize);

    tensor_t result(scalarSize);
    result[nullIndex] = p;
    return result;
  }

  static tensor_t invert(const tensor_t& t)
  {
    return adjugate(t)/determinant(t);
  }

  static expression_t determinant(const tensor_t& t)
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

    return det;
  }
  
  static tensor_t gradToDiv(const tensor_t& operand)
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

  void printTop(const std::string& name) const
  {
    std::ostream& out = std::cout;

    out << name << ":" << std::endl;
    out << valueStack.top();
  }

public:
  FormEvaluationVisitor(const Scenario<dimension>& _scenario, const std::size_t _basisFunctionIndex) : 
    scenario(_scenario), basisFunctionIndex(_basisFunctionIndex)
  {
  }

  expression_t jacobianDeterminant() const
  {
    return determinant(buildJacobian());
  }

  tensor_t getResult() const
  {
    assert(valueStack.size() == 1);
    return valueStack.top();
  }

  tensor_t getTop() const
  {
    assert(!valueStack.empty());
    return valueStack.top();
  }

  virtual void enter(detail::FieldAddition& addition)
  {
    // Nothing to do here
  }

  virtual void exit(detail::FieldAddition& addition)
  {
    const tensor_t second = valueStack.top();
    valueStack.pop();

    const tensor_t first = valueStack.top();
    valueStack.pop();

    valueStack.push(first + second);
    //printTop("FieldAddition");
  }

  virtual void enter(detail::FieldInnerProduct& inner)
  {
    // Nothing to do here
  }

  virtual void exit(detail::FieldInnerProduct& inner)
  {
    const tensor_t second = valueStack.top();
    valueStack.pop();

    const tensor_t first = valueStack.top();
    valueStack.pop();

    valueStack.push(first.inner_product(second));
    //printTop("InnerProduct");
  }

  virtual void enter(detail::FieldOuterProduct& outer)
  {
    // Nothing to do here
  }

  virtual void exit(detail::FieldOuterProduct& outer)
  {
    const tensor_t second = valueStack.top();
    valueStack.pop();

    const tensor_t first = valueStack.top();
    valueStack.pop();

    valueStack.push(first.outer_product(second));
    //printTop("OuterProduct");
  }

  virtual void enter(detail::FieldColonProduct& colon)
  {
    // Nothing to do here
  }

  virtual void exit(detail::FieldColonProduct& colon)
  {
    const tensor_t second = valueStack.top();
    valueStack.pop();

    const tensor_t first = valueStack.top();
    valueStack.pop();

    valueStack.push(first.colon_product(second));
    //printTop("ColonProduct");
  }

  virtual void enter(detail::FieldGradient& gradient)
  {
    // Nothing to do here
  }

  virtual void exit(detail::FieldGradient& gradient)
  {
    const tensor_t value = valueStack.top();
    valueStack.pop();

    valueStack.push(buildGlobalGradient(value));
    //printTop("Gradient");
  }

  virtual void enter(detail::FieldDivergence& divergence)
  {
    // Nothing to do here
  }

  virtual void exit(detail::FieldDivergence& divergence)
  {
    const tensor_t value = valueStack.top();
    valueStack.pop();

    valueStack.push(gradToDiv(buildGlobalGradient(value)));
    //printTop("Divergence");
  }

  virtual void visit(detail::FacetNormal& normal)
  {
    //FIXME: implement me
    CFD_EXCEPTION("Facet normal generation not yet implemented!");
    //printTop("FacetNormal");
  }

  virtual void visit(detail::FieldBasis& basis)
  {
    const FiniteElement<dimension>& element = scenario.getElement(basis.getElement());
    valueStack.push(element.getBasis(basisFunctionIndex, position));
    //printTop("Basis");
  }

  virtual void visit(detail::FieldDiscreteReference& field)
  {
    const FunctionSpaceExpr::expr_ptr functionSpace = field.getDiscreteField().getExpr()->getFunctionSpace();
    const DofMap<dimension>& fieldDofMap = scenario.getDofMap(*functionSpace);

    // This implicitly checks that the DofMap isn't composite.
    const FiniteElement<dimension>& element = *fieldDofMap.getFiniteElement();
    const TensorSize basisSize(element.getRank(), element.getDimension());
    tensor_t fieldValue(basisSize);

    for (std::size_t index=0; index<element.spaceDimension(); ++index)
    {
      fieldValue += element.getBasis(index, position) * 
        ScalarPlaceholder(BasisCoefficient(field.getDiscreteField(), index));
    }
    valueStack.push(fieldValue);
    //printTop("DiscreteField");
  }

  virtual void visit(detail::FieldScalar& s)
  {
    const ScalarPlaceholder scalar(ScalarAccess(s.getValue()));
    valueStack.push(asTensor(scalar));
    //printTop("FieldScalar");
  }
};

}

}

#endif
