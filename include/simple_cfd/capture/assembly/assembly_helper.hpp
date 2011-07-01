#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_ASSEMBLY_HELPER_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_ASSEMBLY_HELPER_HPP

#include <cstddef>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <boost/array.hpp>
#include <simple_cfd/cell_manager.hpp>
#include <simple_cfd/local_assembly_matrix.hpp>
#include <simple_cfd/capture/forms/basis_finder.hpp>
#include "scalar_placeholder.hpp"
#include "position_component.hpp"
#include "form_evaluation_visitor.hpp"

namespace cfd
{

namespace detail
{

namespace
{

template<std::size_t D>
class SymbolicIntegrator
{
private:
  static const std::size_t dimension = D;
  typedef ScalarPlaceholder::expression_t expression_t;
  typedef expression_t::value_map value_map;

  ScalarPlaceholder::expression_t jacobian;
  value_map valueMap;

  static value_map getValueMap(const LocalTransformation<dimension, dimension>& transform)
  {
    value_map valueMap;
    const SmallVector<dimension, expression_t> transformed = transform.getTransformed();

    for(std::size_t d=0; d<dimension; ++d)
      valueMap.bind(ScalarPlaceholder(PositionComponent(d)), transformed[d]);

    return valueMap;
  }

public:
  SymbolicIntegrator(const typename CellManager::ref<dimension>::general cell) : 
    jacobian(cell->getCellReferenceLocalTransformation().getScalingFactor()),
    valueMap(getValueMap(cell->getCellReferenceLocalTransformation()))
  {
  }

  ScalarPlaceholder::expression_t operator()(const ScalarPlaceholder::expression_t& p) const
  {
    expression_t integrand = p.substituteValues(valueMap)*jacobian;
    expression_t::region_t region;

    for(std::size_t d=0; d<dimension; ++d)
      region.setInterval(ScalarPlaceholder(PositionComponent(d)), -1, 1);

    integrand = integrand.integrate(region);
    return integrand.normalised();
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
  typedef ScalarPlaceholder::expression_t expression_t;
  typedef Tensor<dimension, expression_t> tensor_t;

  const Scenario<dimension>& scenario;

public:
  AssemblyHelper(const Scenario<dimension>& _scenario) : scenario(_scenario)
  {
  }

  void assembleBilinearForm(LocalAssemblyMatrix<dimension, expression_t>& matrix, const forms::BilinearForm& form) const
  {
    forms::BasisFinder<dimension> trialFinder(scenario);
    form.getTrialField()->accept(trialFinder);
    const FiniteElement<dimension>* trialElement = trialFinder.getBasis();

    forms::BasisFinder<dimension> testFinder(scenario);
    form.getTestField()->accept(testFinder);
    const FiniteElement<dimension>* testElement = testFinder.getBasis();

    std::vector<tensor_t> trialValues(trialElement->spaceDimension());
    std::vector<tensor_t> testValues(testElement->spaceDimension());

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
        matrix(row, col) += trialValues[trialBasis].colon_product(testValues[testBasis]); 
      }
    }
  }

  LocalAssemblyMatrix<dimension, expression_t> integrate(const LocalAssemblyMatrix<dimension, expression_t>& matrix,
    const MeshEntity& localEntity)
  {
    // TODO: edge integrals, etc.
    const SymbolicIntegrator<dimension> integrator(scenario.getMesh().getReferenceCell());

    // Multiply by local-to-global jacobian
    // TODO: factor jacobian stuff into helper function
    FormEvaluationVisitor<dimension> evaluationVisitor(scenario, 0);
    const expression_t jacobianDet = abs(evaluationVisitor.jacobianDeterminant());
    LocalAssemblyMatrix<dimension, expression_t> result(matrix * jacobianDet);

    // Now change co-ordinates and integrate over reference space, -1 to 1
    std::transform(result.begin(), result.end(), result.begin(), integrator);
    return result;
  }
};

}

}

#endif
