#ifndef SIMPLE_CFD_CAPTURE_SOLVE_OPERATION_HPP
#define SIMPLE_CFD_CAPTURE_SOLVE_OPERATION_HPP

#include <cassert>
#include <map>
#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "fields/fields_fwd.hpp"
#include "fields/named_field.hpp"
#include "fields/field.hpp"
#include "fields/discrete_expr_container.hpp"
#include "fields/discrete_expr_container_builder.hpp"
#include "dimensionless_scenario.hpp"
#include "evaluation/evaluation_strategy.hpp"
#include "exception.hpp"

namespace cfd
{

class SolveOperation
{
private:
  detail::DimensionlessScenario* scenario;
  bool finished;
  std::map<NamedField, Field> newValues;

  //NOTE: we want ownership here, but have to support assignment
  boost::shared_ptr<detail::EvaluationStrategy> evaluationStrategy;

  bool isExecutable() const
  {
    return scenario != NULL && finished;
  }

public:
  SolveOperation() : scenario(NULL), finished(false)
  {
  }

  SolveOperation(detail::DimensionlessScenario& s) : scenario(&s), finished(false)
  {
  }

  void setNewValue(NamedField& f, const Field& value)
  {
    newValues[f] = value;
  }

  void finish()
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
    boost::shared_ptr<detail::EvaluationStrategy> strategy(new detail::EvaluationStrategy(exprContainer, outputExpressions));
    evaluationStrategy.swap(strategy);
  }

  void execute()
  {
    if (!isExecutable())
      CFD_EXCEPTION("Attempted to execute a solve operation before finishing construction of it.");

    scenario->execute(*this);
  }

  template<std::size_t D>
  void executeDimensionTemplated()
  {
    evaluationStrategy->execute<D>();
  }
};

}

#endif
