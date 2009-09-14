#ifndef SIMPLE_CFD_CAPTURE_SOLVE_OPERATION_HPP
#define SIMPLE_CFD_CAPTURE_SOLVE_OPERATION_HPP

#include <cassert>
#include <map>
#include "fields/named_field.hpp"
#include "fields/field.hpp"
#include "fields/discrete_expr_container.hpp"
#include "fields/discrete_expr_container_builder.hpp"

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

    detail::DiscreteExprContainerBuilder containerBuilder;

    for(std::map<NamedField, Field>::const_iterator newValuesIter(newValues.begin()); newValuesIter!=newValues.end(); ++newValuesIter)
    {
      newValuesIter->second.getExpr()->accept(containerBuilder);
    }

    assert(!containerBuilder.containsUndefinedNodes());
    detail::DiscreteExprContainer exprContainer(containerBuilder.getContainer());
  }

  void execute()
  {
    // TODO: Implement me!
  }
};

}

#endif
