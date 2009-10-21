#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_EVALUATION_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_EVALUATION_VISITOR_HPP

#include <cstddef>
#include <cassert>
#include <simple_cfd/capture/fields/discrete_expr_visitor.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class EvaluationVisitor : public DiscreteExprVisitor
{
private:
  static const std::size_t dimension = D;

public:
  // Discrete field related
  virtual void visit(DiscreteFieldElementWise& p)
  {
    assert(false);
  }

  virtual void visit(DiscreteFieldTwoNorm& p)
  {
    assert(false);
  }

  virtual void visit(DiscreteFieldProjection& p)
  {
    assert(false);
  }

  virtual void visit(DiscreteFieldUndefined& u)
  {
    assert(false);
  }

  virtual void visit(DiscreteFieldZero& z)
  {
    assert(false);
  }

  virtual void visit(DiscreteFieldPersistent& p)
  {
    assert(false);
  }


  // Discrete operator related
  virtual void visit(OperatorApplication& a)
  {
    assert(false);
  }

  virtual void visit(OperatorAddition& u)
  {
    assert(false);
  }

  virtual void visit(OperatorAssembly& a)
  {
    assert(false);
  }

  virtual void visit(OperatorUndefined& u)
  {
    assert(false);
  }


  // Scalar related
  virtual void visit(ScalarBinaryOperator& o)
  {
    assert(false);
  }

  virtual void visit(ScalarLiteral& l)
  {
    assert(false);
  }

  virtual void visit(ScalarUndefined& l)
  {
    assert(false);
  }


  // Temporal related
  virtual void visit(DiscreteIndexedScalar& s)
  {
    assert(false);
  }

  virtual void visit(DiscreteIndexedField& s)
  {
    assert(false);
  }

  virtual void visit(DiscreteIndexedOperator& s)
  {
    assert(false);
  }


  // Solve related
  virtual void visit(LinearSolve& s)
  {
    assert(false);
  }
};

}

}

#endif
