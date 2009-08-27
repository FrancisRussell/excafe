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
#include "linear_form.hpp"

namespace cfd
{

namespace forms
{

namespace
{

template<std::size_t D>
class FormEvaluatorVisitor : public FieldVisitor
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
  FormEvaluatorVisitor(const LinearForm& _form, const GeneralCell<dimension>& _cell, const CellVertices<dimension>& _cellVertices,
    const MeshEntity& _localEntity, const vertex<dimension>& v, const Dof<dimension>& _dof) : 
    form(_form), cell(_cell.cloneGeneralCell()), vertices(_cellVertices), localEntity(_localEntity), localVertex(v), 
    dof(_dof), evaluationState(VALUE)
  {
  }

  Tensor<dimension> getResult() const
  {
    assert(valueStack.size() == 1);
    return valueStack.top();
  }

  virtual void enter(FieldAddition& addition)
  {
  }

  virtual void exit(FieldAddition& addition)
  {
    const Tensor<dimension> second = valueStack.top();
    valueStack.pop();

    const Tensor<dimension> first = valueStack.top();
    valueStack.pop();

    valueStack.push(first + second);
  }

  virtual void enter(InnerProduct& inner)
  {
  }

  virtual void exit(InnerProduct& inner)
  {
    const Tensor<dimension> second = valueStack.top();
    valueStack.pop();

    const Tensor<dimension> first = valueStack.top();
    valueStack.pop();

    valueStack.push(first.inner_product(second));
  }

  virtual void enter(OuterProduct& outer)
  {
  }

  virtual void exit(OuterProduct& outer)
  {
    const Tensor<dimension> second = valueStack.top();
    valueStack.pop();

    const Tensor<dimension> first = valueStack.top();
    valueStack.pop();

    valueStack.push(first.outer_product(second));
  }

  virtual void enter(ColonProduct& colon)
  {
  }

  virtual void exit(ColonProduct& colon)
  {
    const Tensor<dimension> second = valueStack.top();
    valueStack.pop();

    const Tensor<dimension> first = valueStack.top();
    valueStack.pop();

    valueStack.push(first.colon_product(second));
  }

  virtual void enter(Gradient& gradient)
  {
    assert(evaluationState == VALUE);
    evaluationState = GRADIENT;
  }

  virtual void exit(Gradient& gradient)
  {
    assert(evaluationState == GRADIENT);
    evaluationState = VALUE;
  }

  virtual void enter(Divergence& divergence)
  {
    assert(evaluationState == VALUE);
    evaluationState = DIVERGENCE;
  }

  virtual void exit(Divergence& divergence)
  {
    assert(evaluationState == DIVERGENCE);
    evaluationState = VALUE;
  }

  virtual void visit(FacetNormal& normal)
  {
    assert(localEntity.getDimension() == dimension-1);
    const Tensor<dimension> facetNormal = cell->getFacetNormal(vertices, localEntity.getIndex(), localVertex);
    valueStack.push(facetNormal);
  }

  virtual void visit(BasisField& basis)
  {
    const FiniteElement<dimension>* const element = 
      boost::any_cast<const FiniteElement<dimension>*>(basis.getElement().getElementPtr());

    assert(element!=NULL);
    valueStack.push(evaluateBasis(*element, dof.getIndex()));
  }

  virtual void visit(DiscreteFieldReference& field)
  {
    const DiscreteField<dimension>* const vector = 
      boost::any_cast<const DiscreteField<dimension>*>(field.getVector().getVectorPtr());

    assert(vector!=NULL);
    assert(!vector->isComposite());

    const FiniteElement<dimension>* const element = vector->getElement();

    std::size_t rank = element->getRank();
    if (evaluationState == DIVERGENCE) --rank;
    if (evaluationState == GRADIENT) ++rank;

    Tensor<dimension> value(rank);

    for(std::size_t i=0; i<element->spaceDimension(); ++i)
    {
      const Dof<dimension> discreteDof(element, dof.getCell(), i);
      double prevTrialCoeff;
      vector->getValues(1, &discreteDof, &prevTrialCoeff);
      value += evaluateBasis(*element, i) * prevTrialCoeff;
    }

    valueStack.push(value);
  }

  virtual void visit(TensorLiteral& literal)
  {
    valueStack.push(boost::any_cast< Tensor<dimension> >(literal.getTensor().getTensor()));
  }
};

}

template<std::size_t D>
class FormEvaluator
{
private:
  static const std::size_t dimension = D;
  std::auto_ptr< GeneralCell<dimension> > cell;
  LinearForm form;

public:
  FormEvaluator(const LinearForm& _form, const GeneralCell<dimension>& _cell) :
    cell(_cell.cloneGeneralCell()), form(_form)
  {
  }

  FormEvaluator(const FormEvaluator& f) : cell(f.cell->cloneGeneralCell()), form(f.form)
  {
  }

  Tensor<dimension> evaluate(const CellVertices<dimension>& vertices, 
    const MeshEntity& localEntity, const vertex<dimension>& v, const Dof<dimension>& dof) const
  {
    FormEvaluatorVisitor<dimension> visitor(form, *cell, vertices, localEntity, v, dof);
    form.accept(visitor);
    return visitor.getResult();
  }
};

}

}

#endif
