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
class PolynomialIntegrator
{
private:
  static const std::size_t dimension = D;
  std::vector< std::pair<vertex<dimension>, double> > rule;

public:
  PolynomialIntegrator(const std::vector<std::pair<vertex<dimension>, double> >& _rule) : rule(_rule)
  {
  }

  ScalarPlaceholder::expression_t operator()(const ScalarPlaceholder::expression_t& p) const
  {
    ScalarPlaceholder::expression_t result;
    ScalarPlaceholder::expression_t::value_map valueMap;

    for(std::size_t point=0; point<rule.size(); ++point)
    {
      const double weight = rule[point].second;

      for(std::size_t d=0; d<dimension; ++d)
      {
        const PositionComponent component(d);
        const ScalarPlaceholder x(component);
        valueMap.bind(x, rule[point].first[d]);
      }
      result += p.substituteValues(valueMap) * weight;
    }
    return result;
  }
};

template<std::size_t D>
class DegreeFinder
{
private:
  static const std::size_t dimension = D;
  boost::array<std::size_t, dimension> degrees;

public:
  DegreeFinder()
  {
    std::fill(degrees.begin(), degrees.end(), 0);
  }

  void operator()(const ScalarPlaceholder::expression_t& p)
  {
    for(std::size_t d=0; d<dimension; ++d)
      degrees[d] = std::max(degrees[d], p.degree(ScalarPlaceholder(PositionComponent(d))));
  }

  boost::array<std::size_t, dimension> getDegrees() const
  {
    return degrees;
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
    DegreeFinder<dimension> degreeFinder;
    degreeFinder = std::for_each(matrix.begin(), matrix.end(), degreeFinder);
    const boost::array<std::size_t, dimension> degrees = degreeFinder.getDegrees();

    const QuadraturePoints<dimension> quadrature(scenario.getMesh().getReferenceCell()->getQuadrature(degrees));

    // FIXME: edge integrals, etc.
    assert(localEntity.getDimension() == dimension);

    const std::vector< std::pair<vertex<dimension>, double> > rule(quadrature.begin(localEntity), quadrature.end(localEntity));
    const PolynomialIntegrator<dimension> integrator(rule);

    // FIXME: factor jacobian stuff into helper function
    FormEvaluationVisitor<dimension> evaluationVisitor(scenario, 0);
    const expression_t jacobianDet = evaluationVisitor.jacobianDeterminant();
    LocalAssemblyMatrix<dimension, expression_t> result(matrix * jacobianDet);
    std::transform(result.begin(), result.end(), result.begin(), integrator);
    return result;
  }
};

}

}

#endif