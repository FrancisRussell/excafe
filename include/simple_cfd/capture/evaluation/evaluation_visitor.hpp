#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_EVALUATION_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_EVALUATION_VISITOR_HPP

#include <cstddef>
#include <cassert>
#include <set>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/capture/capture_fwd.hpp>
#include <simple_cfd/capture/fields/discrete_traits.hpp>
#include <simple_cfd/capture/fields/discrete_field_element_wise.hpp>
#include <simple_cfd/capture/fields/discrete_field_two_norm.hpp>
#include <simple_cfd/capture/fields/discrete_field_zero.hpp>
#include <simple_cfd/capture/fields/discrete_field_projection.hpp>
#include <simple_cfd/capture/fields/discrete_field_apply_bc.hpp>
#include <simple_cfd/capture/fields/operator_assembly.hpp>
#include <simple_cfd/capture/fields/operator_application.hpp>
#include <simple_cfd/capture/fields/operator_addition.hpp>
#include <simple_cfd/capture/fields/operator_apply_bc.hpp>
#include <simple_cfd/capture/fields/scalar_literal.hpp>
#include <simple_cfd/capture/fields/scalar_binary_operator.hpp>
#include <simple_cfd/capture/fields/discrete_indexed_object.hpp>
#include <simple_cfd/capture/fields/linear_solve.hpp>
#include <simple_cfd/discrete_value_traits.hpp>
#include <simple_cfd/numeric/solver.hpp>
#include "discrete_expr_scoping.hpp"
#include "discrete_expr_scoping_visitor.hpp"
#include "expression_values.hpp"

//FIXME: remove me when evil boundary condition hack is over
#include <boost/any.hpp>
#include <simple_cfd/capture/fields/linear_system.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class EvaluationVisitor : public DiscreteExprScopingVisitor
{
private:
  static const std::size_t dimension = D;
  Scenario<dimension>& scenario;
  ExpressionValues<dimension> values;
  ExpressionValues<dimension> kept;

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

  class DiscreteFieldElementWiseEvaluator : public boost::static_visitor<field_value_t>
  {
  private:
    EvaluationVisitor& v;
    DiscreteFieldElementWise& e;

  public:
    DiscreteFieldElementWiseEvaluator(EvaluationVisitor& _v, DiscreteFieldElementWise& _e) : v(_v), e(_e)
    {
    }

    field_value_t operator()(const DiscreteFieldElementWise::add_tag&) const
    {
      return v.getValue(e.getLeft()) + v.getValue(e.getRight());
    }

    field_value_t operator()(const DiscreteFieldElementWise::sub_tag&) const
    {
      return v.getValue(e.getLeft()) - v.getValue(e.getRight());
    }
  };

  class ExecutionHelper : public boost::static_visitor<void>
  {
  private:
    EvaluationVisitor<D>& parent;
    DiscreteExprScoping& currentScope;

  public:
    ExecutionHelper(EvaluationVisitor<D>& _parent, DiscreteExprScoping& _currentScope) : parent(_parent),
      currentScope(_currentScope)
    {
    }

    void operator()(DiscreteExpr* const expr) const
    {
      expr->accept(parent);
    }

    void operator()(TemporalIndexValue* const loopIndex) const
    {
      // FIXME: check sub-scope actually exists
      DiscreteExprScoping& loopScope = currentScope.getLoop(*loopIndex);
      loopScope.accept(parent);
    }
  };

  void execute(DiscreteExprScoping& scope)
  {
    ExecutionHelper helper(*this, scope);
    for(DiscreteExprScoping::const_iterator evalIter(scope.begin()); evalIter!=scope.end(); ++evalIter)
    {
      DiscreteExprScoping::evaluatable_t evaluatable(*evalIter);
      boost::apply_visitor(helper, evaluatable);
    }
  }

  scalar_value_t& getValue(ScalarExpr& e)
  {
    return values.getValue(e);
  }

  field_value_t& getValue(DiscreteFieldExpr& e)
  {
    return values.getValue(e);
  }

  operator_value_t& getValue(OperatorExpr& e)
  {
    return values.getValue(e);
  }

  scalar_value_t getValue(IndexableValue<discrete_scalar_tag>& i, const signed offset)
  {
    if (offset == 0)
    {
      return getValue(*i.getIterationAssignment());
    }
    else if (values.hasValue(i, offset))
    {
      return values.getValue(i, offset);
    }
    else
    {
      return 0;
    }
  }

  field_value_t getValue(IndexableValue<discrete_field_tag>& i, const int offset)
  {
    if (offset == 0)
    {
      return getValue(*i.getIterationAssignment());
    }
    else if (values.hasValue(i, offset))
    {
      return values.getValue(i, offset);
    }
    else
    {
      return DiscreteField<dimension>(scenario.getDofMap(*i.getIterationAssignment()->getFunctionSpace()));
    }
  }

  operator_value_t getValue(IndexableValue<discrete_operator_tag>& i, const int offset)
  {
    if (offset == 0)
    {
      return getValue(*i.getIterationAssignment());
    }
    else if (values.hasValue(i, offset))
    {
      return values.getValue(i, offset);
    }
    else
    {
      return DiscreteOperator<dimension>(scenario.getDofMap(*i.getIterationAssignment()->getTestSpace()),
        scenario.getDofMap(*i.getIterationAssignment()->getTrialSpace()));
    }
  }

  void setValue(ScalarExpr& e, const scalar_value_t& v)
  {
    values.setValue(e, v);
  }

  void setValue(DiscreteFieldExpr& e, const field_value_t& v)
  {
    values.setValue(e, v);
  }

  void setValue(OperatorExpr& e, const operator_value_t& v)
  {
    values.setValue(e, v);
  }

  void setValue(IndexableValue<discrete_scalar_tag>& i, const scalar_value_t& value)
  {
    values.setValue(i, value);
  }

  void setValue(IndexableValue<discrete_field_tag>& i, const field_value_t& value)
  {
    values.setValue(i, value);
  }

  void setValue(IndexableValue<discrete_operator_tag>& i, const operator_value_t& value)
  {
    values.setValue(i, value);
  }

  template<typename discrete_object_tag>
  void resolveIndexableExprsTyped(const std::set<IndexableValue<discrete_object_tag>*>& exprs)
  {

    for(typename std::set<IndexableValue<discrete_object_tag>*>::const_iterator exprIter(exprs.begin()); exprIter!=exprs.end(); 
      ++exprIter)
    {
      IndexableValue<discrete_object_tag>& indexable = **exprIter;
      setValue(indexable, getValue(*indexable.getIterationAssignment()));
    }
  }

  void resolveIndexableExprs(DiscreteExprScoping& scope, TemporalIndexValue* const loopIndex)
  {
    const std::set<IndexableValue<discrete_scalar_tag>*> indexableScalars = loopIndex->getIndexableScalars();
    resolveIndexableExprsTyped(indexableScalars);

    const std::set<IndexableValue<discrete_field_tag>*> indexableFields = loopIndex->getIndexableFields();
    resolveIndexableExprsTyped(indexableFields);

    const std::set<IndexableValue<discrete_operator_tag>*> indexableOperators = loopIndex->getIndexableOperators();
    resolveIndexableExprsTyped(indexableOperators);
  }

public:
  EvaluationVisitor(Scenario<dimension>& _scenario) : 
    scenario(_scenario)
  {
  }

  virtual void visitBlock(DiscreteExprScoping& scope)
  {
    values.enterScope();
    execute(scope);

    if (values.isGlobalScope())
    {
      kept = values;
    }

    values.exitScope();
  }

  virtual void visitLoop(DiscreteExprScoping& scope, TemporalIndexValue* const loopIndex)
  {
    bool looping = true;

    values.enterScope();
    while(looping)
    {
      execute(scope);
      resolveIndexableExprs(scope, loopIndex);

      looping = getValue(loopIndex->getTermination()) == 0;

      if (looping)
        values.completeIteration();
    }
    values.calculateFinals();
    values.exitScope();
  }

  // Discrete field related
  virtual void visit(DiscreteFieldElementWise& p)
  {
    const DiscreteFieldElementWiseEvaluator evaluator(*this, p);
    const DiscreteFieldElementWise::operator_t operation = p.getOperation();
    setValue(p, boost::apply_visitor(evaluator, operation));
  }

  virtual void visit(DiscreteFieldTwoNorm& p)
  {
    const scalar_value_t value = getValue(p.getField()).two_norm();
    setValue(p, value);
  }

  virtual void visit(DiscreteFieldProjection& p)
  {
    const DofMap<dimension> newDofMap(scenario.getDofMap(*p.getFunctionSpace()));
    setValue(p, getValue(p.getField()).project(newDofMap));
  }

  virtual void visit(DiscreteFieldUndefined& u)
  {
    CFD_EXCEPTION("Tried to evaluate an undefined field!");
  }

  virtual void visit(DiscreteFieldZero& z)
  {
    setValue(z, DiscreteField<dimension>(scenario.getDofMap(*z.getFunctionSpace())));
  }

  virtual void visit(DiscreteFieldPersistent& p)
  {
    setValue(p, scenario.getNamedValue(p));
  }

  virtual void visit(DiscreteFieldApplyBC& a)
  {
    //FIXME: implement me!
    assert(false);
    DiscreteField<dimension> newField(getValue(a.getField()));
    //const LinearSystem::load_vector_bc_t bc = boost::any_cast<LinearSystem::load_vector_bc_t>(a.getBoundaryCondition());
    //bc(newField);
    setValue(a, newField);
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
//  const operator_value_t value = getValue(u.getLeft()) + getValue(u.getRight());
//  setValue(u, value);
  }

  virtual void visit(OperatorAssembly& a)
  {
    const DofMap<dimension>& trialDofMap(scenario.getDofMap(*a.getTrialSpace()));
    const DofMap<dimension>& testDofMap(scenario.getDofMap(*a.getTestSpace()));
    DiscreteOperator<dimension> discreteOperator(testDofMap, trialDofMap);
    discreteOperator.assembleForms(scenario, values, a.getBilinearFormIntegralSum());
    setValue(a, discreteOperator);
  }

  virtual void visit(OperatorApplyBC& a)
  {
    //FIXME: implement me!
    assert(false);
    DiscreteOperator<dimension> newOperator(getValue(a.getOperator()));
    //const LinearSystem::system_matrix_bc_t bc = boost::any_cast<LinearSystem::system_matrix_bc_t>(a.getBoundaryCondition());
    //bc(newOperator);
    setValue(a, newOperator);
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
    PETScKrylovSolver solver;
    solver.setMaxIterations(25000);
    solver.setAbsoluteTolerance(1e-4);
    solver.setRelativeTolerance(0.0);
    solver.enablePreconditioner(false);
    
    DiscreteOperator<dimension>& stiffnessMatrix = getValue(s.getOperator());
    DiscreteField<dimension>& loadVector = getValue(s.getField());
    DiscreteField<dimension> unknownVector(stiffnessMatrix.getColMappings());

    solver.solve(stiffnessMatrix.getMatrixHandle(), unknownVector.getVectorHandle(), loadVector.getVectorHandle());

    if (!solver.converged())
    {
      CFD_EXCEPTION("Convergence failure: " + solver.getConvergedReason());
    }
 
    setValue(s, unknownVector);
  }

  ExpressionValues<dimension> getWanted() const
  {
    return kept;
  }
};

}

}

#endif
