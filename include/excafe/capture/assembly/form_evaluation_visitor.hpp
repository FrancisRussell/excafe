#ifndef EXCAFE_CAPTURE_ASSEMBLY_FORM_EVALUATION_VISITOR_HPP
#define EXCAFE_CAPTURE_ASSEMBLY_FORM_EVALUATION_VISITOR_HPP

#include <stack>
#include <algorithm>
#include <excafe/cell_manager.hpp>
#include <excafe/exception.hpp>
#include <excafe/numeric/tensor.hpp>
#include <excafe/numeric/tensor_matrix_view.hpp>
#include <excafe/numeric/polynomial.hpp>
#include <excafe/numeric/functional.hpp>
#include <excafe/capture/forms/linear_form.hpp>
#include <excafe/capture/forms/field_visitor.hpp>
#include <excafe/capture/forms/field_addition.hpp>
#include <excafe/capture/forms/field_inner_product.hpp>
#include <excafe/capture/forms/field_outer_product.hpp>
#include <excafe/capture/forms/field_colon_product.hpp>
#include <excafe/capture/forms/field_gradient.hpp>
#include <excafe/capture/forms/field_divergence.hpp>
#include <excafe/capture/forms/field_basis.hpp>
#include <excafe/capture/forms/field_discrete_reference.hpp>
#include <excafe/capture/forms/field_scalar.hpp>
#include <excafe/capture/forms/facet_normal.hpp>
#include "scalar_placeholder.hpp"
#include "position_placeholder.hpp"
#include "cell_vertices_placeholder.hpp"
#include "scalar_access.hpp"
#include "tensor_operations.hpp"

namespace excafe
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

  TensorOperations<dimension> tensorOps;
  const Scenario<dimension>& scenario;
  PositionPlaceholder position;
  CellVerticesPlaceholder<dimension> cellVertices;
  const std::size_t basisFunctionIndex;
  std::stack<tensor_t> valueStack;

  tensor_t buildGlobalGradient(const tensor_t& operand) const
  {
    const tensor_t localGradient = tensorOps.grad(operand);
    const tensor_t inverseGradient(tensorOps.invert(tensorOps.grad(buildGlobalPosition())));

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
      globalPosition += cell->getCoordinateMapping().getBasis(i, position) * cellVertices[i];

    return globalPosition;
  }

  tensor_t buildJacobian() const
  {
    return tensorOps.transpose(tensorOps.grad(buildGlobalPosition()));
  }

  static tensor_t asTensor(const expression_t& p)
  {
    const TensorSize scalarSize(0, dimension);
    const TensorIndex nullIndex(scalarSize);

    tensor_t result(scalarSize);
    result[nullIndex] = p;
    return result;
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
    return tensorOps.determinant(buildJacobian());
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

    valueStack.push(tensorOps.gradToDiv(buildGlobalGradient(value)));
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
