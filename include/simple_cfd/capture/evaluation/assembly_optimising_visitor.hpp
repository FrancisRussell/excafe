#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_ASSEMBLY_OPTIMISING_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_ASSEMBLY_OPTIMISING_VISITOR_HPP

#include <cstddef>
#include <set>
#include <iostream>
#include <boost/foreach.hpp>
#include <simple_cfd/finite_element.hpp>
#include <simple_cfd/dof_map.hpp>
#include <simple_cfd/local_assembly_matrix.hpp>
#include <simple_cfd/numeric/functional.hpp>
#include <simple_cfd/capture/fields/operator_assembly.hpp>
#include <simple_cfd/capture/fields/discrete_expr_visitor.hpp>
#include <simple_cfd/capture/forms/bilinear_form_integral_sum.hpp>
#include <simple_cfd/capture/assembly/scalar_placeholder.hpp>
#include <simple_cfd/capture/assembly/assembly_helper.hpp>
#include <simple_cfd/codegen/ufc_evaluator.hpp>
#include <simple_cfd/capture/evaluation/local_assembly_matrix_evaluator.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class AssemblyOptimisingVisitor : public DiscreteExprVisitor
{
private:
  static const std::size_t dimension = D;
  Scenario<dimension>& scenario;

public:
  AssemblyOptimisingVisitor(Scenario<dimension>& _scenario) : scenario(_scenario)
  {
  }
  
  // Discrete field related
  virtual void visit(DiscreteFieldElementWise& p) {}
  virtual void visit(DiscreteFieldTwoNorm& p) {}
  virtual void visit(DiscreteFieldProjection& p) {}
  virtual void visit(DiscreteFieldUndefined& u) {}
  virtual void visit(DiscreteFieldZero& z) {}
  virtual void visit(DiscreteFieldPersistent& p) {}
  virtual void visit(DiscreteFieldApplyBC& a) {}

  // Discrete operator related
  virtual void visit(OperatorApplication& a) {}
  virtual void visit(OperatorAddition& u) {}
  virtual void visit(OperatorUndefined& u) {}
  virtual void visit(OperatorApplyBC& a) {}

  // Scalar related
  virtual void visit(ScalarBinaryOperator& o) {}
  virtual void visit(ScalarLiteral& l) {}
  virtual void visit(ScalarUndefined& l) {}

  // Temporal related
  virtual void visit(DiscreteIndexedScalar& s) {}
  virtual void visit(DiscreteIndexedField& s) {}
  virtual void visit(DiscreteIndexedOperator& s) {}

  // Solve related
  virtual void visit(LinearSolve& s) {}

  virtual void visit(OperatorAssembly& a)
  {
    typedef ScalarPlaceholder::expression_t expression_t;
    typedef LocalAssemblyMatrix<dimension, expression_t> local_matrix_t;

    const forms::BilinearFormIntegralSum sum = a.getBilinearFormIntegralSum();
    const local_matrix_t localMatrix = scenario.constructCellIntegralAssemblyMatrix(sum);
    const LocalAssemblyMatrixEvaluator<dimension> evaluator = 
      codegen::UFCEvaluator<dimension>::construct(scenario, localMatrix);

    /* 
       FIXME: We only save a cell integral. We need to support saving multiple optimised integrals
       depending on whether they're cell or facet, and if facet, internal or external facets.
    */
    scenario.setOptimisedCellIntegral(a, evaluator);
  }
};

}

}

#endif
