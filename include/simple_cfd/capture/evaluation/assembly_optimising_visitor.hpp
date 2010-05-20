#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_ASSEMBLY_OPTIMISING_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_ASSEMBLY_OPTIMISING_VISITOR_HPP

#include <cstddef>
#include <set>
#include <iostream>
#include <simple_cfd/finite_element.hpp>
#include <simple_cfd/dof_map.hpp>
#include <simple_cfd/local_assembly_matrix.hpp>
#include <simple_cfd/numeric/functional.hpp>
#include <simple_cfd/capture/fields/operator_assembly.hpp>
#include <simple_cfd/capture/fields/discrete_expr_visitor.hpp>
#include <simple_cfd/capture/forms/bilinear_form_integral_sum.hpp>
#include <simple_cfd/capture/assembly/scalar_placeholder.hpp>
#include <simple_cfd/capture/assembly/assembly_helper.hpp>

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
    typedef expression_t::optimised_t optimised_expression_t;
    typedef LocalAssemblyMatrix<dimension, expression_t> local_matrix_t;
    typedef LocalAssemblyMatrix<dimension, optimised_expression_t> opt_local_matrix_t;
    typedef FiniteElement<dimension> finite_element_t;

    const DofMap<dimension>& testMapping = scenario.getDofMap(*a.getTestSpace());
    const DofMap<dimension>& trialMapping = scenario.getDofMap(*a.getTrialSpace());

    const std::set<const finite_element_t*> trialElements(trialMapping.getFiniteElements());
    const std::set<const finite_element_t*> testElements(testMapping.getFiniteElements());
    const forms::BilinearFormIntegralSum sum = a.getBilinearFormIntegralSum();

    //FIXME: Hard coded to cell integrals
    const forms::BilinearFormIntegralSum::const_iterator sumBegin = sum.begin_dx();
    const forms::BilinearFormIntegralSum::const_iterator sumEnd = sum.end_dx();

    AssemblyHelper<dimension> assemblyHelper(scenario);
    local_matrix_t localMatrix(testElements, trialElements);
    
    for(forms::BilinearFormIntegralSum::const_iterator formIter = sumBegin; formIter!=sumEnd; ++formIter)
    {
      assemblyHelper.assembleBilinearForm(localMatrix, *formIter);
      std::cout << "Assembled local-matrix expression " << 1 + formIter - sumBegin << " of " << sumEnd - sumBegin << std::endl;
    }

    //FIXME: Hard coded to cell integrals
    const MeshEntity localCellEntity(dimension, 0);
    localMatrix = assemblyHelper.integrate(localMatrix, localCellEntity);

    std::cout << "Integrated local-matrix expression..." << std::endl;

    //std::cout << localMatrix << std::endl;

    const opt_local_matrix_t optimisedLocalMatrix(localMatrix.transform(PolynomialOptimiser<expression_t>()));

    std::cout << "Built optimised local-matrix expression..." << std::endl;

    /* 
       FIXME: We only save a cell integral. We need to support saving multiple optimised integrals
       depending on whether they're cell or facet, and if facet, internal or external facets.
    */
    scenario.setOptimisedCellIntegral(a, optimisedLocalMatrix);
  }
};

}

}

#endif
