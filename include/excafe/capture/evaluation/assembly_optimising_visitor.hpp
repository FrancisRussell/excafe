#ifndef EXCAFE_CAPTURE_EVALUATION_ASSEMBLY_OPTIMISING_VISITOR_HPP
#define EXCAFE_CAPTURE_EVALUATION_ASSEMBLY_OPTIMISING_VISITOR_HPP

#include <cstddef>
#include <set>
#include <iostream>
#include <boost/foreach.hpp>
#include <excafe/finite_element.hpp>
#include <excafe/dof_map.hpp>
#include <excafe/local_assembly_matrix.hpp>
#include <excafe/numeric/functional.hpp>
#include <excafe/capture/fields/operator_assembly.hpp>
#include <excafe/capture/fields/discrete_expr_visitor.hpp>
#include <excafe/capture/forms/bilinear_form_integral_sum.hpp>
#include <excafe/capture/assembly/scalar_placeholder.hpp>
#include <excafe/capture/assembly/assembly_helper.hpp>
#include <excafe/codegen/ufc_evaluator.hpp>
#include <excafe/capture/evaluation/local_assembly_matrix_evaluator.hpp>
#include <excafe/capture/evaluation/local_assembly_matrix_interpreter.hpp>

namespace excafe
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

    const bool interpret = false;
    if (interpret)
    {
      const LocalAssemblyMatrixEvaluator<dimension> evaluator =
        LocalAssemblyMatrixInterpreter<dimension>::construct(scenario, localMatrix);

      scenario.setOptimisedCellIntegral(a, evaluator);
    }
    else
    {
      const LocalAssemblyMatrixEvaluator<dimension> evaluator = 
        codegen::UFCEvaluator<dimension>::construct(scenario, localMatrix);

      /* 
         FIXME: We only save a cell integral. We need to support saving multiple optimised integrals
         depending on whether they're cell or facet, and if facet, internal or external facets.
      */
      scenario.setOptimisedCellIntegral(a, evaluator);
    }
  }
};

}

}

#endif
