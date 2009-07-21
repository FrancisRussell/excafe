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
  const CellVertices<dimension> vertices;
  const vertex<dimension> localVertex;
  const Dof<dimension> dof;


  EvaluationState evaluationState;
  std::stack< Tensor<dimension> > valueStack;

public:
  FormEvaluatorVisitor(const LinearForm& _form, const CellVertices<dimension>& _cellVertices,
    const vertex<dimension>& v, const Dof<dimension>& _dof) : form(_form), vertices(_cellVertices),
    localVertex(v), dof(_dof), evaluationState(VALUE)
  {
  }

  Tensor<dimension> getResult() const
  {
    assert(valueStack.size() == 1);
    return valueStack.top();
  }

  virtual void enter(Addition& addition)
  {
  }

  virtual void exit(Addition& addition)
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

  virtual void visit(BasisField& basis)
  {
    const FiniteElement<dimension>* const element = 
      boost::any_cast<const FiniteElement<dimension>*>(basis.getElement().getElementPtr());

  }

  virtual void visit(DiscreteField& field) = 0;
  virtual void visit(TensorLiteral& literal) = 0;
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
    const vertex<dimension>& v, const Dof<dimension>& dof)
  {
    FormEvaluatorVisitor<dimension> visitor(vertices, v, dof);
    form.accept(visitor);
    return visitor.getResult();
  }
};

}

}

#endif
