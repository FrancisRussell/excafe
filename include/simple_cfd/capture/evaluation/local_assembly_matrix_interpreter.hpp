#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_LOCAL_ASSEMBLY_MATRIX_INTERPRETER_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_LOCAL_ASSEMBLY_MATRIX_INTERPRETER_HPP

#include <cstddef>
#include <algorithm>
#include <boost/utility.hpp>
#include <boost/foreach.hpp>
#include <simple_cfd/numeric/functional.hpp>
#include "local_assembly_matrix_evaluator_impl.hpp"
#include "local_assembly_matrix_evaluator.hpp"

namespace cfd
{

namespace detail
{

template<std::size_t D>
class LocalAssemblyMatrixInterpreter : public LocalAssemblyMatrixEvaluatorImpl<D>,
                                       boost::noncopyable
{
public:
  static const std::size_t dimension = D;
  typedef ScalarPlaceholder::expression_t expression_t;

private:
  const Scenario<dimension>& scenario;
  const LocalAssemblyMatrix<dimension, expression_t> localAssemblyMatrix;
  std::set<ScalarPlaceholder> placeholders;

  LocalAssemblyMatrixInterpreter(const Scenario<dimension>& _scenario,
    const LocalAssemblyMatrix<dimension, expression_t>& _localAssemblyMatrix) : 
    scenario(_scenario), localAssemblyMatrix(_localAssemblyMatrix)
  {
    ExpressionVariableCollector<expression_t> collector;
    collector = std::for_each(localAssemblyMatrix.begin(), localAssemblyMatrix.end(), collector);
    placeholders = collector.getVariables();
  }

public:
  static LocalAssemblyMatrixEvaluator<dimension> construct(const Scenario<dimension>& scenario,
    const LocalAssemblyMatrix<dimension, expression_t>& localAssemblyMatrix)
  {
    std::auto_ptr< LocalAssemblyMatrixEvaluatorImpl<dimension> > 
      interpreter(new LocalAssemblyMatrixInterpreter(scenario, localAssemblyMatrix));

    return LocalAssemblyMatrixEvaluator<dimension>(interpreter);
  }

  void evaluate(LocalAssemblyMatrix<dimension, double>& matrix,
                std::size_t cid,
                const ExpressionValues<dimension>& values) const
  {
    const ScalarPlaceholderEvaluator<dimension> placeholderEvaluator(scenario, values, cid);
    expression_t::value_map valueMap;

    BOOST_FOREACH(const ScalarPlaceholder& placeholder, placeholders)
    {
      valueMap.bind(placeholder, placeholderEvaluator(placeholder));
    }

    const ExpressionEvaluator<expression_t> evaluator(valueMap);
    std::transform(localAssemblyMatrix.begin(), localAssemblyMatrix.end(), matrix.begin(), evaluator);
  }
};

}

}

#endif
