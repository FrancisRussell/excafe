#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_EVALUATION_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_EVALUATION_VISITOR_HPP

#include <cstddef>
#include <cassert>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <simple_cfd/exception.hpp>
//#include <simple_cfd/capture/scenario.hpp>
#include <simple_cfd/capture/capture_fwd.hpp>
#include <simple_cfd/capture/fields/discrete_expr_visitor.hpp>
#include <simple_cfd/capture/fields/discrete_traits.hpp>
#include <simple_cfd/capture/fields/discrete_field_element_wise.hpp>
#include <simple_cfd/capture/fields/discrete_field_two_norm.hpp>
#include <simple_cfd/capture/fields/discrete_field_zero.hpp>
#include <simple_cfd/capture/fields/operator_application.hpp>
#include <simple_cfd/capture/fields/operator_addition.hpp>
#include <simple_cfd/capture/fields/scalar_literal.hpp>
#include <simple_cfd/capture/fields/scalar_binary_operator.hpp>
#include <simple_cfd/capture/fields/discrete_indexed_object.hpp>
#include <simple_cfd/discrete_value_traits.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class EvaluationVisitor : public DiscreteExprVisitor
{
private:
  static const std::size_t dimension = D;
  Scenario<dimension>& scenario;

  typedef typename DiscreteValueTraits<discrete_scalar_tag, D>::value_t scalar_value_t;
  typedef typename DiscreteValueTraits<discrete_field_tag, D>::value_t field_value_t;
  typedef typename DiscreteValueTraits<discrete_operator_tag, D>::value_t operator_value_t;
  
  class ScalarBinaryOperatorEvaluator : public boost::static_visitor<scalar_value_t>
  {
  private:
    const scalar_value_t left;
    const scalar_value_t right;

  public:
    ScalarBinaryOperatorEvaluator(const scalar_value_t _left, const scalar_value_t _right) : left(_left), right(_right)
    {
    }

    scalar_value_t operator()(const ScalarBinaryOperator::add_tag&) const { return left+right; }
    scalar_value_t operator()(const ScalarBinaryOperator::sub_tag&) const { return left-right; }
    scalar_value_t operator()(const ScalarBinaryOperator::div_tag&) const { return left/right; }
    scalar_value_t operator()(const ScalarBinaryOperator::mul_tag&) const { return left*right; }
    scalar_value_t operator()(const ScalarBinaryOperator::lt_tag&) const  { return left<right ? 1.0 : 0.0; }
    scalar_value_t operator()(const ScalarBinaryOperator::gt_tag&) const  { return left>right ? 1.0 : 0.0; }
    scalar_value_t operator()(const ScalarBinaryOperator::lte_tag&) const { return left<=right ? 1.0 : 0.0; }
    scalar_value_t operator()(const ScalarBinaryOperator::gte_tag&) const { return left>=right ? 1.0 : 0.0; }
    scalar_value_t operator()(const ScalarBinaryOperator::eq_tag&) const  { return left==right ? 1.0 : 0.0; }
  };

  scalar_value_t& getValue(ScalarExpr& e)
  {
  }

  field_value_t& getValue(DiscreteFieldExpr& e)
  {
  }

  operator_value_t& getValue(OperatorExpr& e)
  {
  }

  scalar_value_t& getValue(IndexableValue<discrete_scalar_tag>& i, const signed offset)
  {
  }

  field_value_t& getValue(IndexableValue<discrete_field_tag>& i, const int offset)
  {
  }

  operator_value_t& getValue(IndexableValue<discrete_operator_tag>& i, const int offset)
  {
  }

  void setValue(ScalarExpr& e, const scalar_value_t& v)
  {
  }

  void setValue(DiscreteFieldExpr& e, const field_value_t& v)
  {
  }

  void setValue(OperatorExpr& e, const operator_value_t& v)
  {
  }

public:
  EvaluationVisitor(Scenario<dimension>& _scenario) : scenario(_scenario)
  {
  }

  // Discrete field related
  virtual void visit(DiscreteFieldElementWise& p)
  {
    const field_value_t value = getValue(p.getLeft()) + getValue(p.getRight());
    setValue(p, value);
  }

  virtual void visit(DiscreteFieldTwoNorm& p)
  {
    const scalar_value_t value = getValue(p.getField()).two_norm();
    setValue(p, value);
  }

  virtual void visit(DiscreteFieldProjection& p)
  {
    assert(false);
  }

  virtual void visit(DiscreteFieldUndefined& u)
  {
    CFD_EXCEPTION("Tried to evaluate an undefined field!");
  }

  virtual void visit(DiscreteFieldZero& z)
  {
    setValue(z, DiscreteField<dimension>(scenario.getDofMap(*z.getFunctionSpace())));
    assert(false);
  }

  virtual void visit(DiscreteFieldPersistent& p)
  {
    assert(false);
  }


  // Discrete operator related
  virtual void visit(OperatorApplication& a)
  {
    const field_value_t value = getValue(a.getOperator()) * getValue(a.getField());
    setValue(a, value);
  }

  virtual void visit(OperatorAddition& u)
  {
    assert(false && "We don't know how to add operators!");
//    const operator_value_t value = getValue(u.getLeft()) + getValue(u.getRight());
//    setValue(u, value);
  }

  virtual void visit(OperatorAssembly& a)
  {
    assert(false);
  }

  virtual void visit(OperatorUndefined& u)
  {
    CFD_EXCEPTION("Tried to evaluate an undefined operator!");
  }


  // Scalar related
  virtual void visit(ScalarBinaryOperator& o)
  {
    const scalar_value_t left = getValue(o.getLeft());
    const scalar_value_t right = getValue(o.getRight());
    const ScalarBinaryOperatorEvaluator evaluator(left, right);
    const ScalarBinaryOperator::operator_t operation = o.getOperator();
    const scalar_value_t value = boost::apply_visitor(evaluator, operation);
    setValue(o, value);
  }

  virtual void visit(ScalarLiteral& l)
  {
    const scalar_value_t value = l.getValue();
    setValue(l, value);
  }

  virtual void visit(ScalarUndefined& l)
  {
    CFD_EXCEPTION("Tried to evaluate an undefined scalar!");
  }


  // Temporal related
  virtual void visit(DiscreteIndexedScalar& s)
  {
    if (s.isInsideLoop())
    {
      setValue(s, getValue(s.getParent(), s.getOffsetValue()));
    }
  }

  virtual void visit(DiscreteIndexedField& s)
  {
    if (s.isInsideLoop())
    {
      setValue(s, getValue(s.getParent(), s.getOffsetValue()));
    }
  }

  virtual void visit(DiscreteIndexedOperator& s)
  {
    if (s.isInsideLoop())
    {
      setValue(s, getValue(s.getParent(), s.getOffsetValue()));
    }
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
