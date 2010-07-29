#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_DUAL_FORM_EVALUATOR_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_DUAL_FORM_EVALUATOR_HPP

#include <string>
#include <cassert>
#include <sstream>
#include <map>
#include <utility>
#include <iosfwd>
#include <boost/foreach.hpp>
#include <simple_cfd/dof.hpp>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/cell_vertices.hpp>
#include <simple_cfd/cell_manager.hpp>
#include <simple_cfd/mesh_entity.hpp>
#include <simple_cfd/numeric/functional.hpp>
#include <simple_cfd/numeric/tensor.hpp>
#include <simple_cfd/capture/forms/form_evaluator_visitor.hpp>
#include <simple_cfd/capture/assembly/scalar_placeholder.hpp>
#include <simple_cfd/capture/assembly/form_evaluation_visitor.hpp>
#include <simple_cfd/capture/assembly/scalar_placeholder_evaluator.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class DualFormEvaluator : public FieldVisitor
{
private:
  static const std::size_t dimension = D;
  typedef typename CellManager::ref<dimension>::general cell_ref_t;
  const cell_ref_t cell;
  const CellVertices<dimension> vertices;
  const MeshEntity localEntity;
  const vertex<dimension> localVertex;
  const Dof<dimension> dof;
  const Scenario<dimension>& scenario;
  const detail::ExpressionValues<dimension>& values;
  
  enum EvaluationState
  {
    VALUE,
    GRADIENT,
    DIVERGENCE
  };

  static bool tensorsEqual(const Tensor<dimension, double>& a, const Tensor<dimension, double>& b) 
  {
    Tensor<dimension, double> tensorDelta(a);
    tensorDelta-=b;

    double abs = 0.0;
    BOOST_FOREACH(const double s, tensorDelta)
    {
      abs += std::abs(s);
    }

    return abs < 1e-3;
  }

  cfd::forms::FormEvaluatorVisitor<dimension> concreteVisitor;
  cfd::detail::FormEvaluationVisitor<dimension> captureVisitor;
  EvaluationState evaluationState;
  
  // TODO: remove duplicated code from DiscreteOperator
  std::map<ScalarPlaceholder, double> evaluatePlaceholders(const std::set<detail::ScalarPlaceholder>& placeholders) const
  {
    const ScalarPlaceholderEvaluator<dimension> evaluator(scenario, values, dof.getCell(), localVertex);
    std::map<ScalarPlaceholder, double> values;

    BOOST_FOREACH(const ScalarPlaceholder& placeholder, placeholders)
    {
      values.insert(std::make_pair(placeholder, evaluator(placeholder)));
    }

    return values;
  }

  void compareTops()
  {
    if (evaluationState == VALUE)
    {
      const Tensor<dimension, double> concreteTop = concreteVisitor.getTop();
      const Tensor<dimension, assembly_polynomial_t> captureTop = captureVisitor.getTop();

      if (concreteTop.getSize() != captureTop.getSize())
      {
        std::ostringstream error;
        error << "Tensor size mismatch when comparing visitor states." << std::endl;
        error << "Concrete visitor had tensor size of " << concreteTop.getSize() << "." << std::endl;
        error << "Capturing visitor had tensor size of " << captureTop.getSize() << "." << std::endl;
        CFD_EXCEPTION(error.str());
      }

      typedef assembly_polynomial_t::optimised_t optimised_t;
      const Tensor<dimension, optimised_t> optimisedTensor = captureTop.transform(PolynomialOptimiser<assembly_polynomial_t>());
      
      PolynomialVariableCollector<optimised_t> collector;
      collector = std::for_each(optimisedTensor.begin(), optimisedTensor.end(), collector);
      const std::map<ScalarPlaceholder, double> placeholderValues(evaluatePlaceholders(collector.getVariables()));

      const PolynomialEvaluator<optimised_assembly_polynomial_t> evaluator(placeholderValues);
      const Tensor<dimension, double> evaluatedTensor = optimisedTensor.transform(evaluator);

      if (!tensorsEqual(concreteTop, evaluatedTensor))
      {
        std::cout << "Comparing visitor states... " << std::endl;
        std::cout << "Position: " << localVertex << std::endl;
        std::cout << "DoF: " << dof << std::endl;
        std::cout << "Capture visitor: " << captureTop << std::endl;
        std::cout << "Concrete visitor: " << concreteTop << std::endl;
        std::cout << "Evaluated Capture visitor: " << evaluatedTensor << std::endl;

        CFD_EXCEPTION("Exiting after tensor mismatch.");
      }
    }
    else
    {
      std::cout << "Cannot compare visitor states due to differing div & grad implementations." << std::endl;
    }
  }

  template<typename T>
  void dualEntry(const std::string& name, T& t)
  {
    concreteVisitor.enter(t);
    captureVisitor.enter(t);
    std::cout << "Operation: " << name << " entry." << std::endl;
  }
  
  template<typename T>
  void dualExit(const std::string& name, T& t)
  {
    concreteVisitor.exit(t);
    captureVisitor.exit(t);
    std::cout << "Operation: " << name << " exit." << std::endl;
    compareTops();
  }

  template<typename T>
  void dualVisit(const std::string& name, T& t)
  {
    concreteVisitor.visit(t);
    captureVisitor.visit(t);
    std::cout << "Operation: " << name << std::endl;
    compareTops();
  }

public:
  DualFormEvaluator(const cell_ref_t _cell, const Scenario<dimension>& _scenario,
    const detail::ExpressionValues<dimension>& _values,
    const CellVertices<dimension>& _cellVertices, const MeshEntity& _localEntity,
    const vertex<dimension>& v, const Dof<dimension>& _dof) :
    cell(_cell), vertices(_cellVertices), localEntity(_localEntity), localVertex(v),
    dof(_dof), scenario(_scenario), values(_values),
    concreteVisitor(cell, scenario, values, vertices, localEntity, localVertex, dof),
    captureVisitor(scenario, dof.getIndex()), evaluationState(VALUE)
  {
  }

  Tensor<dimension, double> getResult() const
  {
    return concreteVisitor.getResult();
  }

  // Non Terminals
  virtual void enter(FieldAddition& addition)
  {
    dualEntry("FieldAddition", addition);
  }

  virtual void exit(FieldAddition& addition)
  {
    dualExit("FieldAddition", addition);
  }

  virtual void enter(FieldInnerProduct& inner) 
  {
    dualEntry("FieldInnerProduct", inner);
  }

  virtual void exit(FieldInnerProduct& inner)
  {
    dualExit("FieldInnerProduct", inner);
  }

  virtual void enter(FieldOuterProduct& outer)
  {
    dualEntry("FieldOuterProduct", outer);
  }

  virtual void exit(FieldOuterProduct& outer)
  {
    dualExit("FieldOuterProduct", outer);
  }

  virtual void enter(FieldColonProduct& colon)
  {
    dualEntry("FieldColonProduct", colon);
  }

  virtual void exit(FieldColonProduct& colon)
  {
    dualExit("FieldColonProduct", colon);
  }

  virtual void enter(FieldGradient& gradient)
  {
    assert(evaluationState == VALUE);
    evaluationState = GRADIENT;

    dualEntry("FieldGradient", gradient);
  }

  virtual void exit(FieldGradient& gradient)
  {
    assert(evaluationState == GRADIENT);
    evaluationState = VALUE;

    dualExit("FieldGradient", gradient);
  }

  virtual void enter(FieldDivergence& divergence)
  {
    assert(evaluationState == VALUE);
    evaluationState = DIVERGENCE;

    dualEntry("FieldDivergence", divergence);
  }

  virtual void exit(FieldDivergence& divergence)
  {
    assert(evaluationState == DIVERGENCE);
    evaluationState = VALUE;

    dualExit("FieldDivergence", divergence);
  }

  // Terminals
  virtual void visit(FacetNormal& normal)
  {
    dualVisit("FacetNormal", normal);
  }

  virtual void visit(FieldBasis& basis)
  {
    dualVisit("FieldBasis", basis);
  }

  virtual void visit(FieldDiscreteReference& field)
  {
    dualVisit("FieldDiscreteReference", field);
  }

  virtual void visit(FieldScalar& s)
  {
    dualVisit("FieldScalar", s);
  }
};

}

}

#endif
