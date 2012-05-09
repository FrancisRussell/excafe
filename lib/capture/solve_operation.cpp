#include <cassert>
#include <boost/shared_ptr.hpp>
#include <excafe/capture/solve_operation.hpp>
#include <excafe/capture/dimensionless_scenario.hpp>
#include <excafe/capture/fields/discrete_expr_container.hpp>
#include <excafe/capture/fields/discrete_expr_container_builder.hpp>
#include <excafe/capture/fields/discrete_expr_container_builder.hpp>
#include <excafe/capture/evaluation/assembly_optimising_visitor.hpp>
#include <excafe/exception.hpp>

namespace excafe
{

bool SolveOperation::isExecutable() const
{
  return scenario != NULL && finished;
}

SolveOperation::SolveOperation() : scenario(NULL), finished(false)
{
}

SolveOperation::SolveOperation(detail::DimensionlessScenario& s) : scenario(&s), finished(false)
{
}

void SolveOperation::setNewValue(NamedField& f, const Field& value)
{
  newValues[f] = value;
}

void SolveOperation::finish()
{
  finished=true;

  std::set<detail::DiscreteExpr*> outputExpressions;
  detail::DiscreteExprContainerBuilder containerBuilder;

  for(std::map<NamedField, Field>::const_iterator newValuesIter(newValues.begin()); newValuesIter!=newValues.end(); ++newValuesIter)
  {
    newValuesIter->second.getExpr()->accept(containerBuilder);
    outputExpressions.insert(&(*newValuesIter->second.getExpr()));
  }

  assert(!containerBuilder.containsUndefinedNodes());
  detail::DiscreteExprContainer exprContainer(containerBuilder.getContainer());

  scenario->resolveFunctionSpaces(exprContainer.getFunctionSpaces());
  scenario->optimiseLocalAssemblies(exprContainer.getOperatorExpressions());
  boost::shared_ptr<detail::EvaluationStrategy> strategy(new detail::EvaluationStrategy(exprContainer, outputExpressions));
  evaluationStrategy.swap(strategy);
}

void SolveOperation::execute()
{
  if (!isExecutable())
    CFD_EXCEPTION("Attempted to execute a solve operation before finishing construction of it.");

   scenario->execute(*this);
}


}
