#ifndef SIMPLE_CFD_FORMS_FORM_EVAULATOR_HPP
#define SIMPLE_CFD_FORMS_FORM_EVALUATOR_HPP

#include <cstddef>
#include <stack>
#include <cassert>
#include <memory>
#include <boost/any.hpp>
#include <simple_cfd/numeric/tensor.hpp>
#include <simple_cfd/cell_vertices.hpp>
#include <simple_cfd/vertex.hpp>
#include <simple_cfd/general_cell.hpp>
#include <simple_cfd/mesh_entity.hpp>
#include <simple_cfd/capture/capture_fwd.hpp>
#include <simple_cfd/capture/evaluation/evaluation_fwd.hpp>
#include "linear_form.hpp"
#include "field_visitor.hpp"
#include "field_addition.hpp"
#include "field_inner_product.hpp"
#include "field_outer_product.hpp"
#include "field_colon_product.hpp"
#include "field_gradient.hpp"
#include "field_divergence.hpp"
#include "field_basis.hpp"
#include "field_discrete_reference.hpp"
#include "field_scalar.hpp"
#include "facet_normal.hpp"

namespace cfd
{

namespace forms
{

namespace
{

template<std::size_t D>
class FormEvaluatorVisitor : public detail::FieldVisitor
{
private:
  enum EvaluationState
  {
    VALUE,
    GRADIENT,
    DIVERGENCE
  };

  static const std::size_t dimension = D;
  const LinearForm form;
  const std::auto_ptr< GeneralCell<dimension> > cell;
  const CellVertices<dimension> vertices;
  const MeshEntity localEntity;
  const vertex<dimension> localVertex;
  const Dof<dimension> dof;
  Scenario<dimension>& scenario;
  detail::ExpressionValues<dimension>& values;


  EvaluationState evaluationState;
  std::stack< Tensor<dimension> > valueStack;

  Tensor<dimension> evaluateBasis(const FiniteElement<dimension>& element, const std::size_t index) const
  {
    if (evaluationState == GRADIENT)
    {
      return element.evaluateGradient(vertices, index, localVertex);
    }
    else if (evaluationState == DIVERGENCE)
    {
      return element.evaluateDivergence(vertices, index, localVertex);
    }
    else
    {
      return element.evaluateTensor(vertices, index, localVertex);
    }
  }

public:
  FormEvaluatorVisitor(const GeneralCell<dimension>& _cell, Scenario<dimension>& _scenario,
    detail::ExpressionValues<dimension>& _values, const LinearForm& _form, 
    const CellVertices<dimension>& _cellVertices, const MeshEntity& _localEntity, 
    const vertex<dimension>& v, const Dof<dimension>& _dof) : 
    form(_form), cell(_cell.cloneGeneralCell()), vertices(_cellVertices), localEntity(_localEntity), localVertex(v), 
    dof(_dof), scenario(_scenario), values(_values), evaluationState(VALUE)
  {
  }

  Tensor<dimension> getResult() const
  {
    assert(valueStack.size() == 1);
    return valueStack.top();
  }

  virtual void enter(detail::FieldAddition& addition)
  {
  }

  virtual void exit(detail::FieldAddition& addition)
  {
    const Tensor<dimension> second = valueStack.top();
    valueStack.pop();

    const Tensor<dimension> first = valueStack.top();
    valueStack.pop();

    valueStack.push(first + second);
  }

  virtual void enter(detail::FieldInnerProduct& inner)
  {
  }

  virtual void exit(detail::FieldInnerProduct& inner)
  {
    const Tensor<dimension> second = valueStack.top();
    valueStack.pop();

    const Tensor<dimension> first = valueStack.top();
    valueStack.pop();

    valueStack.push(first.inner_product(second));
  }

  virtual void enter(detail::FieldOuterProduct& outer)
  {
  }

  virtual void exit(detail::FieldOuterProduct& outer)
  {
    const Tensor<dimension> second = valueStack.top();
    valueStack.pop();

    const Tensor<dimension> first = valueStack.top();
    valueStack.pop();

    valueStack.push(first.outer_product(second));
  }

  virtual void enter(detail::FieldColonProduct& colon)
  {
  }

  virtual void exit(detail::FieldColonProduct& colon)
  {
    const Tensor<dimension> second = valueStack.top();
    valueStack.pop();

    const Tensor<dimension> first = valueStack.top();
    valueStack.pop();

    valueStack.push(first.colon_product(second));
  }

  virtual void enter(detail::FieldGradient& gradient)
  {
    assert(evaluationState == VALUE);
    evaluationState = GRADIENT;
  }

  virtual void exit(detail::FieldGradient& gradient)
  {
    assert(evaluationState == GRADIENT);
    evaluationState = VALUE;
  }

  virtual void enter(detail::FieldDivergence& divergence)
  {
    assert(evaluationState == VALUE);
    evaluationState = DIVERGENCE;
  }

  virtual void exit(detail::FieldDivergence& divergence)
  {
    assert(evaluationState == DIVERGENCE);
    evaluationState = VALUE;
  }

  virtual void visit(detail::FacetNormal& normal)
  {
    assert(localEntity.getDimension() == dimension-1);
    const Tensor<dimension> facetNormal = cell->getFacetNormal(vertices, localEntity.getIndex(), localVertex);
    valueStack.push(facetNormal);
  }

  virtual void visit(detail::FieldBasis& basis)
  {
    const FiniteElement<dimension>& element = scenario.getElement(basis.getElement());
    valueStack.push(evaluateBasis(element, dof.getIndex()));
  }

  virtual void visit(detail::FieldDiscreteReference& field)
  {

    const DiscreteField<dimension>& vector = values.getValue(*field.getDiscreteField().getExpr());
    assert(!vector.isComposite());

    const FiniteElement<dimension>* const element = vector.getElement();

    std::size_t rank = element->getRank();
    if (evaluationState == DIVERGENCE) --rank;
    if (evaluationState == GRADIENT) ++rank;

    Tensor<dimension> value(rank);

    for(std::size_t i=0; i<element->spaceDimension(); ++i)
    {
      const Dof<dimension> discreteDof(element, dof.getCell(), i);
      double prevTrialCoeff;
      vector.getValues(1, &discreteDof, &prevTrialCoeff);
      value += evaluateBasis(*element, i) * prevTrialCoeff;
    }

    valueStack.push(value);
  }

  virtual void visit(detail::FieldScalar& s)
  {
    const std::size_t rank = 0;
    Tensor<dimension> value(rank);
    value[NULL] = values.getValue(*s.getValue().getExpr());
    valueStack.push(value);
  }

/*
  TODO: Where the hell did this go?

  virtual void visit(TensorLiteral& literal)
  {
    valueStack.push(boost::any_cast< Tensor<dimension> >(literal.getTensor().getTensor()));
  }
*/
};

}

template<std::size_t D>
class FormEvaluator
{
private:
  static const std::size_t dimension = D;
  std::auto_ptr< GeneralCell<dimension> > cell;
  Scenario<dimension>* scenario;
  detail::ExpressionValues<dimension>* values;
  LinearForm form;

public:
  FormEvaluator(const GeneralCell<dimension>& _cell, Scenario<dimension>& _scenario, 
    detail::ExpressionValues<dimension>& _values, const LinearForm& _form) :
    cell(_cell.cloneGeneralCell()), scenario(&_scenario), values(&_values), form(_form)
  {
  }

  FormEvaluator(const FormEvaluator& f) : cell(f.cell->cloneGeneralCell()), scenario(f.scenario),
    values(f.values), form(f.form)
  {
  }

  Tensor<dimension> evaluate(const CellVertices<dimension>& vertices, 
    const MeshEntity& localEntity, const vertex<dimension>& v, const Dof<dimension>& dof) const
  {
    FormEvaluatorVisitor<dimension> visitor(*cell, *scenario, *values, form, vertices, localEntity, v, dof);
    form.accept(visitor);
    return visitor.getResult();
  }
};

}

}

#endif
