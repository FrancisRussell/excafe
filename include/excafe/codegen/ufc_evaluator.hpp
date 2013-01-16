#ifndef EXCAFE_CODEGEN_UFC_EVALUATOR_HPP
#define EXCAFE_CODEGEN_UFC_EVALUATOR_HPP

#include <map>
#include <set>
#include <memory>
#include <sstream>
#include <iomanip>
#include <ufc.h>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <excafe/cse/cse_optimiser.hpp>
#include <excafe/exception.hpp>
#include <excafe/capture/capture_fwd.hpp>
#include <excafe/capture/fields/function_space_expr.hpp>
#include <excafe/local_assembly_matrix.hpp>
#include <excafe/capture/assembly/scalar_placeholder.hpp>
#include <excafe/capture/assembly/scalar_placeholder_evaluator.hpp>
#include <excafe/capture/evaluation/local_assembly_matrix_evaluator.hpp>
#include <excafe/capture/evaluation/local_assembly_matrix_evaluator_impl.hpp>
#include <excafe/util/timer.hpp>
#include "ufc_integral_generator.hpp"
#include "dynamic_cxx.hpp"

namespace excafe
{

namespace codegen
{

namespace detail
{

class FieldCollector : public boost::static_visitor<void>
{
private:
  std::set<Field::expr_ptr> fields;
  std::set<Scalar::expr_ptr> scalars;

public:
  void operator()(const excafe::detail::PositionComponent& c)
  {
  }

  void operator()(const excafe::detail::CellVertexComponent& c)
  {
  }

  void operator()(const excafe::detail::ScalarAccess& s)
  {
    scalars.insert(s.getExpr());
  }

  void operator()(const excafe::detail::BasisCoefficient& c)
  {
    fields.insert(c.getField());
  }

  void operator()(const excafe::detail::GenericSymbol& s)
  {
    CFD_EXCEPTION("GenericSymbol found in ScalarPlaceholder. This should never happen.");
  }

  void operator()(const boost::blank&)
  {
    CFD_EXCEPTION("boost::blank found in ScalarPlaceholder. This should never happen.");
  }

  std::set<Field::expr_ptr> getFields() const
  {
    return fields;
  }

  std::set<Scalar::expr_ptr> getScalars() const
  {
    return scalars;
  }
};

struct UFCContext
{
  boost::scoped_array<double> cellVertexValues;
  boost::scoped_array<double*> cellVertexPointers;

  boost::scoped_array<double> coefficientValues;
  boost::scoped_array<double*> coefficientPointers;
};

}

template<std::size_t D>
class UFCEvaluator : public excafe::detail::LocalAssemblyMatrixEvaluatorImpl<D>,
                     public boost::noncopyable
{
public:
  static const std::size_t dimension = D;
  typedef excafe::detail::ScalarPlaceholder::expression_t expression_t;

private:
  const Scenario<dimension>& scenario;
  const excafe::detail::LocalAssemblyMatrix<dimension, expression_t>& localAssemblyMatrix;
  std::string className;

  // References to the DSO and dynamically constructed ufc::cell_integral.
  boost::scoped_ptr<DynamicCXX> dsoHandler;
  boost::scoped_ptr<ufc::cell_integral> cellIntegral;

  // The number of entries in each UFC coefficent.
  std::vector<std::size_t> coefficientSizes;

  // The ScalarPlaceholders that need to be evaluated to form the UFC coefficient array.
  std::vector<excafe::detail::ScalarPlaceholder> coefficientPlaceholders;

  // These hold the ordering of fields and scalars as passed via the UFC interface.
  UFCIntegralGenerator::coefficient_index_map_t coefficientIndices;

  // Struct containing information that will be passed to tabulate_tensor
  detail::UFCContext context;

  void addCoefficient(const boost::variant<Field::expr_ptr, Scalar::expr_ptr>& coefficient,
                      const std::vector<excafe::detail::ScalarPlaceholder>& placeholders)
  {
    coefficientIndices.insert(std::make_pair(coefficient, coefficientIndices.size()));
    coefficientSizes.push_back(placeholders.size());
    coefficientPlaceholders.insert(coefficientPlaceholders.end(), placeholders.begin(), placeholders.end());
  }

  void addCoefficientField(const Field::expr_ptr& fieldExpr)
  {
    using excafe::detail::ScalarPlaceholder;
    using excafe::detail::BasisCoefficient;

    const excafe::detail::FunctionSpaceExpr::expr_ptr functionSpaceExpr = fieldExpr->getFunctionSpace();
    const DofMap<dimension>& dofMap = scenario.getDofMap(*functionSpaceExpr);
    const std::size_t dofs = dofMap.getFiniteElement()->spaceDimension();

    std::vector<ScalarPlaceholder> placeholders;
    for(std::size_t dof=0; dof<dofs; ++dof)
      placeholders.push_back(ScalarPlaceholder(BasisCoefficient(fieldExpr, dof)));

    addCoefficient(fieldExpr, placeholders);
  }

  void addCoefficientScalar(const Scalar::expr_ptr& scalarExpr)
  {
    using excafe::detail::ScalarPlaceholder;

    const Scalar scalar(scalarExpr);
    const excafe::detail::ScalarAccess scalarAccess(scalar);
    const std::vector<ScalarPlaceholder> placeholders(1, ScalarPlaceholder(scalarAccess));

    addCoefficient(scalarExpr, placeholders);
  }

  UFCEvaluator(const Scenario<dimension>& _scenario,
               const excafe::detail::LocalAssemblyMatrix<dimension, expression_t>& _localAssemblyMatrix) : 
    scenario(_scenario), localAssemblyMatrix(_localAssemblyMatrix)
  {
  }

  void computeOrdering()
  {
    using excafe::detail::ScalarPlaceholder;

    // Collect all unknown symbols in local assembly matrix.
    ExpressionVariableCollector<expression_t> collector;
    BOOST_FOREACH(const expression_t& e, localAssemblyMatrix)
      collector(e);

    const std::set<ScalarPlaceholder> unknowns = collector.getVariables();

    // Isolate references to fields and scalars.
    detail::FieldCollector fieldCollector;
    BOOST_FOREACH(const ScalarPlaceholder& unknown, unknowns)
      unknown.apply(fieldCollector);

    const std::set<Field::expr_ptr> fields = fieldCollector.getFields();
    const std::set<Scalar::expr_ptr> scalars = fieldCollector.getScalars();

    // For now, we arrange fields in an arbitrary order, followed by scalars.
    BOOST_FOREACH(const Field::expr_ptr& field, fields)
      addCoefficientField(field);
      
    BOOST_FOREACH(const Scalar::expr_ptr& scalar, scalars)
      addCoefficientScalar(scalar);
  }

  void setClassName(const std::string name)
  {
    className = name;
  }

  std::string getInstantiatorName() const
  {
    return std::string("new") + className;
  }

  std::string getCode()
  {
    using excafe::detail::ScalarPlaceholder;

    std::vector<expression_t> expressions;

    unsigned index = 0;
    BOOST_FOREACH(const expression_t& expr, localAssemblyMatrix)
    {
      std::cout << "Adding local assembly matrix polynomial: " << index++ << std::endl;
      const expression_t normalised(expr.normalised());
      expressions.push_back(normalised);
    }

    std::cout << "Calling CSE..." << std::endl;

    util::Timer timer;
    timer.start();
    cse::CSEOptimiser<ScalarPlaceholder> optimiser(expressions.begin(), expressions.end());
    timer.stop();

    std::ostringstream source;
    source << "#include <cassert>\n";
    source << "#include <cmath>\n";
    source << "#include <ufc.h>\n\n";

    const long minutes = timer.getIntegerSeconds()/60;
    const double seconds = static_cast<double>(timer.getIntegerSeconds()%60) + timer.getFractionalSeconds();

    source << "// Common sub-expression elimination pass took ";
    source << minutes << " minute" << (minutes != 1 ? "s":"") << " and ";
    source << std::fixed << std::setprecision(2);
    source << seconds << " seconds (wall clock).\n\n";

    UFCIntegralGenerator generator(source, coefficientIndices);
    generator.outputPrefix();
    optimiser.accept(generator);
    generator.outputPostfix();

    // Generate method used to instantiate new class.
    setClassName(generator.getClassName());
    source << "\n";
    source << "extern \"C\" ufc::cell_integral* " << getInstantiatorName() << "()\n";
    source << "{\n";
    source << "  return new " << className << "();\n";
    source << "}\n";

    return source.str();
  }

  void buildCodeAndOpen()
  {
    const std::string code = getCode();
    dsoHandler.reset(new DynamicCXX(code));
    dsoHandler->compileAndLoad();

    typedef ufc::cell_integral* (*instantiator_t)();
    const instantiator_t instantiator =
      reinterpret_cast<instantiator_t>(dsoHandler->getFunction(getInstantiatorName()));

    cellIntegral.reset(instantiator());
  }

  void initialiseContext()
  {
    const Mesh<dimension>& mesh = scenario.getMesh();
    const std::size_t verticesPerCell = mesh.getReferenceCell()->numEntities(0);

    context.cellVertexValues.reset(new double[dimension*verticesPerCell]);
    context.cellVertexPointers.reset(new double*[verticesPerCell]);

    context.coefficientValues.reset(new double[coefficientPlaceholders.size()]);
    context.coefficientPointers.reset(new double*[coefficientSizes.size()]);

    double* vertexPtr = context.cellVertexValues.get();
    for(std::size_t v=0; v<verticesPerCell; ++v)
    {
      context.cellVertexPointers[v] = vertexPtr;
      vertexPtr += dimension;
    }

    double* coefficientPtr = context.coefficientValues.get();
    for(std::size_t i=0; i<coefficientSizes.size(); ++i)
    {
      context.coefficientPointers[i] = coefficientPtr;
      coefficientPtr += coefficientSizes[i];
    }
  }

  void evaluate(excafe::detail::LocalAssemblyMatrix<dimension, double>& matrix,
                const std::size_t cid,
                const excafe::detail::ExpressionValues<dimension>& values) const
  {
    const CellVertices<dimension> vertices(scenario.getMesh().getCoordinates(cid));
    for(std::size_t v=0; v<vertices.size(); ++v)
    {
      for(std::size_t d=0; d<dimension; ++d)
      {
        context.cellVertexValues[v*dimension + d] = vertices[v][d];
      }
    }

    const excafe::detail::ScalarPlaceholderEvaluator<dimension> evaluator(scenario, values, cid);
    std::transform(coefficientPlaceholders.begin(), 
                   coefficientPlaceholders.end(),
                   context.coefficientValues.get(), 
                   evaluator);

    ufc::cell cell;
    cell.topological_dimension = cell.geometric_dimension = dimension;
    cell.coordinates = context.cellVertexPointers.get();

    cellIntegral->tabulate_tensor(matrix.data(), context.coefficientPointers.get(), cell);
  }

public:
  static excafe::detail::LocalAssemblyMatrixEvaluator<dimension> construct(const Scenario<dimension>& scenario, 
    const excafe::detail::LocalAssemblyMatrix<dimension, expression_t>& localAssemblyMatrix)
  {
    std::auto_ptr<UFCEvaluator> evaluator(new UFCEvaluator(scenario, localAssemblyMatrix));
    evaluator->computeOrdering();
    evaluator->buildCodeAndOpen();
    evaluator->initialiseContext();

    std::auto_ptr< excafe::detail::LocalAssemblyMatrixEvaluatorImpl<dimension> > impl(evaluator);
    return excafe::detail::LocalAssemblyMatrixEvaluator<dimension>(impl);
  }

  static std::string getCode(const Scenario<dimension>& scenario, 
    const excafe::detail::LocalAssemblyMatrix<dimension, expression_t>& localAssemblyMatrix)
  {
    std::auto_ptr<UFCEvaluator> evaluator(new UFCEvaluator(scenario, localAssemblyMatrix));
    evaluator->computeOrdering();
    return evaluator->getCode();
  }

  ~UFCEvaluator()
  {
    // We must delete the run-time generated ufc::cell_integral before closing the DSO.
    cellIntegral.reset();
    dsoHandler.reset();
  }
};

}

}

#endif
