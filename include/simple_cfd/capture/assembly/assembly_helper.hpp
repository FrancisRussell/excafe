#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_ASSEMBLY_HELPER_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_ASSEMBLY_HELPER_HPP

#include <cstddef>
#include <vector>
#include <utility>
#include <algorithm>
#include <simple_cfd/cell_manager.hpp>
#include <simple_cfd/local_assembly_matrix.hpp>
#include <simple_cfd/capture/forms/basis_finder.hpp>
#include "scalar_placeholder.hpp"
#include "position_component.hpp"
#include "form_evaluation_visitor.hpp"
#include "assembly_polynomial.hpp"

namespace cfd
{

namespace detail
{

namespace
{

template<std::size_t D>
class PolynomialIntegrator
{
private:
  static const std::size_t dimension = D;
  std::vector< std::pair<vertex<dimension>, double> > rule;

public:
  PolynomialIntegrator(const std::vector<std::pair<vertex<dimension>, double> >& _rule) : rule(_rule)
  {
  }

  assembly_polynomial_t operator()(const assembly_polynomial_t& p) const
  {
    assembly_polynomial_t result;
    for(std::size_t point=0; point<rule.size(); ++point)
    {
      const double weight = rule[point].second;
      assembly_polynomial_t weighted = p*weight;

      for(std::size_t d=0; d<dimension; ++d)
      {
        const PositionComponent component(d);
        const ScalarPlaceholder dx(component);
        weighted = weighted.substituteValue(dx, rule[point].first[d]);
      }
      result += weight;
    }
    return result;
  }
};

}

template<std::size_t D>
class AssemblyHelper
{
public:
  static const std::size_t dimension = D;

private:
  typedef typename CellManager::ref<dimension>::general cell_ref_t;
  typedef assembly_polynomial_t polynomial_t;

  const Scenario<dimension>& scenario;

public:
  AssemblyHelper(const Scenario<dimension>& _scenario) : scenario(_scenario)
  {
  }

  void assembleBilinearForm(LocalAssemblyMatrix<dimension, polynomial_t>& matrix, const forms::BilinearForm& form) const
  {
    forms::BasisFinder<dimension> trialFinder(scenario);
    form.getTrialField()->accept(trialFinder);
    const FiniteElement<dimension>* trialElement = trialFinder.getBasis();

    forms::BasisFinder<dimension> testFinder(scenario);
    form.getTestField()->accept(testFinder);
    const FiniteElement<dimension>* testElement = testFinder.getBasis();

    std::vector<polynomial_t> trialValues(trialElement->spaceDimension());
    std::vector<polynomial_t> testValues(testElement->spaceDimension());

    for(std::size_t trialBasis=0; trialBasis < trialValues.size(); ++trialBasis)
    {
      FormEvaluationVisitor<dimension> evaluationVisitor(scenario, trialBasis);
      form.getTrialField()->accept(evaluationVisitor);
      trialValues[trialBasis] = evaluationVisitor.getResult();
    }

    for(std::size_t testBasis=0; testBasis < testValues.size(); ++testBasis)
    {
      FormEvaluationVisitor<dimension> evaluationVisitor(scenario, testBasis);
      form.getTestField()->accept(evaluationVisitor);
      testValues[testBasis] = evaluationVisitor.getResult();
    }
    
    for(std::size_t trialBasis=0; trialBasis < trialValues.size(); ++trialBasis)
    {
      for(std::size_t testBasis=0; testBasis < testValues.size(); ++testBasis)
      {
        const std::size_t row = matrix.getTestOffset(*testElement, testBasis);
        const std::size_t col = matrix.getTrialOffset(*trialElement, trialBasis);

        matrix(row, col) += trialValues[trialBasis] * testValues[testBasis]; 
      }
    }
  }

  LocalAssemblyMatrix<dimension, polynomial_t> integrate(const LocalAssemblyMatrix<dimension, polynomial_t>& matrix,
    const MeshEntity& localEntity)
  {
    const QuadraturePoints<dimension> quadrature(scenario.getMesh().getReferenceCell()->getQuadrature(5));

    // FIXME: edge integrals, etc.
    assert(localEntity.getDimension() == dimension);

    const std::vector< std::pair<vertex<dimension>, double> > rule(quadrature.begin(localEntity), quadrature.end(localEntity));
    const PolynomialIntegrator<dimension> integrator(rule);

    // FIXME: factor jacobian stuff into helper function
    FormEvaluationVisitor<dimension> evaluationVisitor(scenario, 0);
    const polynomial_t jacobianDet = evaluationVisitor.jacobianDeterminant();
    LocalAssemblyMatrix<dimension, polynomial_t> result(matrix * jacobianDet);
    std::transform(result.begin(), result.end(), result.begin(), integrator);
    return result;
  }
};

}

}

#endif
