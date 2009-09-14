#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_CONTAINER_BUILDER_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_CONTAINER_BUILDER_HPP

#include "discrete_expr_container.hpp"
#include "discrete_expr_visitor.hpp"
#include "discrete_field_element_wise.hpp"
#include "discrete_field_two_norm.hpp"
#include "discrete_field_projection.hpp"
#include "discrete_field_undefined.hpp"
#include "discrete_field_zero.hpp"
#include "discrete_field_persistent.hpp"
#include "operator_application.hpp"
#include "operator_addition.hpp"
#include "operator_assembly.hpp"
#include "operator_undefined.hpp"
#include "scalar_binary_operator.hpp"
#include "scalar_literal.hpp"
#include "scalar_undefined.hpp"
#include "linear_solve.hpp"
#include "discrete_indexed_object.hpp"
#include "indexable_value.hpp"

namespace cfd
{

namespace detail
{

class DiscreteExprContainerBuilderFormVisitor : public FieldVisitor
{
public:

};

class DiscreteExprContainerBuilder : public DiscreteExprVisitor
{
private:
  bool undefinedNodes;
  DiscreteExprContainer container;

  template<typename discrete_object_tag>
  void handleIndexedNode(IndexableValue<discrete_object_tag>& indexableValue)
  {
    if (container.insert(indexableValue))
    {
      handleTemporalIndex(*indexableValue.getIndexVariable());
      indexableValue.getIterationAssignment()->accept(*this);

      for(typename IndexableValue<discrete_object_tag>::init_iterator initIter(indexableValue.begin_inits());
          initIter!=indexableValue.end_inits(); ++initIter)
      {
        initIter->accept(*this);
      }
    }
  }

  void handleTemporalIndex(TemporalIndexValue& v)
  {
    if (container.insert(v))
    {
      v.getTermination().accept(*this);
    }
  }

public:
  DiscreteExprContainerBuilder() : undefinedNodes(false)
  {
  }

  // Discrete field related
  virtual void visit(DiscreteFieldElementWise& p)
  {
    if (container.insert(p))
    {
      p.getLeft().accept(*this);
      p.getRight().accept(*this);
    }
  }

  virtual void visit(DiscreteFieldTwoNorm& p)
  {
    if (container.insert(p))
    {
      p.getField().accept(*this);
    }
  }

  virtual void visit(DiscreteFieldProjection& p)
  {
    if (container.insert(p))
    {
      p.getField().accept(*this);
    }
  }

  virtual void visit(DiscreteFieldUndefined& u)
  {
    undefinedNodes = true;
    container.insert(u);
  }

  virtual void visit(DiscreteFieldZero& z)
  {
    container.insert(z);
  }

  virtual void visit(DiscreteFieldPersistent& p)
  {
    container.insert(p);
  }

  // Discrete operator related
  virtual void visit(OperatorApplication& a)
  {
    if (container.insert(a))
    {
      a.getOperator().accept(*this);
      a.getField().accept(*this);
    }
  }

  virtual void visit(OperatorAddition& u)
  {
    if (container.insert(u))
    {
      u.getLeft().accept(*this);
      u.getRight().accept(*this);
    }
  }

  virtual void visit(OperatorAssembly& a)
  {
    // TODO: handle referenced fields in the assembled expressions
    container.insert(a);
  }

  virtual void visit(OperatorUndefined& u)
  {
    undefinedNodes = true;
    container.insert(u);
  }

  // Scalar related
  virtual void visit(ScalarBinaryOperator& o)
  {
    if (container.insert(o))
    {
      o.getLeft().accept(*this);
      o.getRight().accept(*this);
    }
  }

  virtual void visit(ScalarLiteral& l)
  {
    container.insert(l);
  }

  virtual void visit(ScalarUndefined& u)
  {
    undefinedNodes = true;
    container.insert(u);
  }

  // Temporal related
  virtual void visit(DiscreteIndexedScalar& s)
  {
    handleIndexedNode(s.getParent());
  }

  virtual void visit(DiscreteIndexedField& s)
  {
    handleIndexedNode(s.getParent());
  }

  virtual void visit(DiscreteIndexedOperator& s)
  {
    handleIndexedNode(s.getParent());
  }

  // Solve related
  virtual void visit(LinearSolve& s)
  {
    if (container.insert(s))
    {
      s.getOperator().accept(*this);
      s.getField().accept(*this);
    }
  }
};

}

}

#endif
