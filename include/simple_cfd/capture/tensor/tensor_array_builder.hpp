#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_BUILDER_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_BUILDER_HPP

#include <cstddef>
#include <stack>
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

  TensorArrayHelper<dimension> helper;
  std::stack<TensorArrayRef> valueStack;

  TensorArrayRef pop()
  {
    const TensorArrayRef top = valueStack.top();
    valueStack.pop();
    return top;
  }

public:
  TensorArrayBuilder(const cell_ref_t _cell) : helper(_cell)
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

  virtual void exit(FieldGradient& gradient)
  {
    //FIXME: implement me!
    assert(false);
  }

  virtual void enter(FieldDivergence& divergence)
  {
  }

  virtual void exit(FieldDivergence& divergence)
  {
    //FIXME: implement me!
    assert(false);
  }

  // Terminals
  virtual void visit(FacetNormal& normal)
  {
    //FIXME: implement me!
    assert(false);
  }

  virtual void visit(FieldBasis& basis)
  {
    //FIXME: implement me!
    assert(false);
  }

  virtual void visit(FieldDiscreteReference& field)
  {
    //FIXME: implement me!
    assert(false);
  }

  virtual void visit(FieldScalar& s)
  {
    //FIXME: implement me!
    assert(false);
  }
};

}

}

#endif
