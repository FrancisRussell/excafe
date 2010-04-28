#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_EVALUATOR_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_EVALUATOR_HPP

#include <cstddef>
#include <cassert>
#include "scalar_placeholder.hpp"
#include <simple_cfd/exception.hpp>
#include <simple_cfd/dof.hpp>
#include <simple_cfd/finite_element.hpp>
#include <simple_cfd/capture/evaluation/expression_values.hpp>

namespace cfd
{

namespace detail
{

namespace
{

/* 
  We need to stick the operator()s in this class, otherwise the 
  operator()(const ScalarPlaceholder&) matches them all
*/
template<std::size_t D>
class ScalarPlaceholderEvaluatorHelper
{
private:
  static const std::size_t dimension = D;
  const Scenario<dimension>& scenario;
  const ExpressionValues<dimension>& values;
  const std::size_t cid;

public:
  typedef double result_type;

  ScalarPlaceholderEvaluatorHelper(const Scenario<dimension>& _scenario, const ExpressionValues<dimension>& _values, 
    const std::size_t _cid) : scenario(_scenario), values(_values), cid(_cid)
  {
  }

  result_type operator()(const PositionComponent& c) const
  {
    CFD_EXCEPTION("ScalarPlaceholderEvaluator currently doesn't handle position placeholders.");
  }

  result_type operator()(const CellVertexComponent& c) const
  {
    const Mesh<dimension>& mesh(scenario.getMesh());
    const CellVertices<dimension> vertices(mesh.getCoordinates(cid));
    return vertices[c.getVertexID()][c.getComponent()];
  }

  result_type operator()(const ScalarAccess& c) const
  {
    return values.getValue(*c.getExpr());
  }

  result_type operator()(const BasisCoefficient& c) const
  {
    const DiscreteField<dimension>& field = values.getValue(*c.getField());
    assert(!field.isComposite());

    const FiniteElement<dimension>* const element = field.getElement();
    const Dof<dimension> discreteDof(element, cid, c.getIndex());

    double coeff;
    field.getValues(1, &discreteDof, &coeff);
    return coeff;
  }
};

}

template<std::size_t D>
class ScalarPlaceholderEvaluator
{
private:
  static const std::size_t dimension = D;
  const ScalarPlaceholderEvaluatorHelper<dimension> helper;

public:
  ScalarPlaceholderEvaluator(const Scenario<dimension>& _scenario, const ExpressionValues<dimension>& _values,
    const std::size_t _cid) : helper(_scenario, _values, _cid)
  {
  }

  double operator()(const ScalarPlaceholder& s) const
  {
    return s.apply(helper);
  }
};

}

}

#endif
