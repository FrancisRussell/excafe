#ifndef SIMPLE_CFD_CAPTURE_SOLVE_OPERATION_HPP
#define SIMPLE_CFD_CAPTURE_SOLVE_OPERATION_HPP

#include <cassert>
#include <map>
#include "fields/fields_fwd.hpp"
#include "fields/named_field.hpp"
#include "fields/field.hpp"
#include "fields/discrete_expr_container.hpp"
#include "fields/discrete_expr_container_builder.hpp"
#include "evaluation/evaluation_strategy.hpp"

namespace cfd
{

class SolveOperation
{
private:
  std::map<NamedField, Field> newValues;

public:
  SolveOperation()
  {
  }

  void setNewValue(NamedField& f, const Field& value)
  {
    newValues[f] = value;
  }

  void finish()
  {
    // TODO: Implement me!

    std::set<detail::DiscreteExpr*> outputExpressions;
    detail::DiscreteExprContainerBuilder containerBuilder;

    for(std::map<NamedField, Field>::const_iterator newValuesIter(newValues.begin()); newValuesIter!=newValues.end(); ++newValuesIter)
    {
      newValuesIter->second.getExpr()->accept(containerBuilder);
      outputExpressions.insert(&(*newValuesIter->second.getExpr()));
    }

    assert(!containerBuilder.containsUndefinedNodes());
    detail::DiscreteExprContainer exprContainer(containerBuilder.getContainer());
    detail::EvaluationStrategy evaluationStrategy(exprContainer, outputExpressions);
  }

  void execute()
  {
    // TODO: Implement me!
  }
};

}

#endif
