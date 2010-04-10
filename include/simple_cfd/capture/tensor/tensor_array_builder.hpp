#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_BUILDER_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_BUILDER_HPP

#include <cstddef>
#include <stack>
#include "tensor_placeholder.hpp"
#include "tensor_array_helper.hpp"
#include <simple_cfd/capture/forms/field_visitor.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class TensorArrayBuilder : public FieldVisitor
{
private:
  static const std::size_t dimension = D;
  typedef typename CellManager::ref<dimension>::general cell_ref_t;

  TensorArrayHelper helper;
  cell_ref_t cell;
  TensorPlaceholder position;
  TensorPlaceholder cellVertices;
  std::stack<TensorArrayRef> valueStack;

  TensorArrayRef pop()
  {
    const TensorArrayRef top = valueStack.top();
    valueStack.pop();
    return top;
  }

public:
  TensorArrayBuilder(const cell_ref_t _cell, const TensorPlaceholder& _position,
    const TensorPlaceholder& _cellVertices) :
    cell(_cell), position(_position), cellVertices(_cellVertices)
  {
  }

  virtual void enter(FieldAddition& addition)
  {
  }

  virtual void exit(FieldAddition& addition)
  {
    const TensorArrayRef second = pop();
    const TensorArrayRef first = pop();
    valueStack.push(helper.add(first, second));
  }

  virtual void enter(FieldInnerProduct& inner)
  {
  }

  virtual void exit(FieldInnerProduct& inner)
  {
    const TensorArrayRef second = pop();
    const TensorArrayRef first = pop();
    valueStack.push(helper.inner(first, second));
  }

  virtual void enter(FieldOuterProduct& outer)
  {
  }

  virtual void exit(FieldOuterProduct& outer)
  {
    const TensorArrayRef second = pop();
    const TensorArrayRef first = pop();
    valueStack.push(helper.outer(first, second));
  }

  virtual void enter(FieldColonProduct& colon)
  {
  }

  virtual void exit(FieldColonProduct& colon)
  {
    const TensorArrayRef second = pop();
    const TensorArrayRef first = pop();
    valueStack.push(helper.colon(first, second));
  }

  virtual void enter(FieldGradient& gradient)
  {
  }

  virtual void exit(FieldGradient& gradient) = 0;

  virtual void enter(FieldDivergence& divergence)
  {
  }

  virtual void exit(FieldDivergence& divergence) = 0;

  // Terminals
  virtual void visit(FacetNormal& normal) = 0;
  virtual void visit(FieldBasis& basis) = 0;
  virtual void visit(FieldDiscreteReference& field) = 0;
  virtual void visit(FieldScalar& s) = 0;
};

}

}

#endif
