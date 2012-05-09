#ifndef EXCAFE_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_EVALUATOR_HPP
#define EXCAFE_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_EVALUATOR_HPP

#include <cstddef>
#include <cassert>
#include <boost/optional.hpp>
#include "scalar_placeholder.hpp"
#include <excafe/exception.hpp>
#include <excafe/vertex.hpp>
#include <excafe/dof.hpp>
#include <excafe/finite_element.hpp>
#include <excafe/capture/evaluation/expression_values.hpp>

namespace excafe
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
  const boost::optional< vertex<dimension> > position;

public:
  typedef double result_type;

  ScalarPlaceholderEvaluatorHelper(const Scenario<dimension>& _scenario, 
    const ExpressionValues<dimension>& _values, 
    const std::size_t _cid, 
    const boost::optional< vertex<dimension> >& _position) : 
    scenario(_scenario), values(_values), cid(_cid), position(_position)
  {
  }

  result_type operator()(const PositionComponent& c) const
  {
    if (!position)
    {
      CFD_EXCEPTION("ScalarPlaceholderEvaluator hasn't been supplied a local co-ordinate.");
    }
    else
    {
      return (*position)[c.getComponent()];
    }
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

  result_type operator()(const GenericSymbol& s) const
  {
    CFD_EXCEPTION("GenericSymbol found in ScalarPlaceholder. This should never happen.");
  }

  result_type operator()(const boost::blank&) const
  {
    CFD_EXCEPTION("boost::blank found in ScalarPlaceholder. This should never happen.");
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
  ScalarPlaceholderEvaluator(const Scenario<dimension>& _scenario, 
    const ExpressionValues<dimension>& _values,
    const std::size_t _cid, 
    const boost::optional< vertex<dimension> >& _position = boost::optional< vertex<dimension> >()) : 
    helper(_scenario, _values, _cid, _position)
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
