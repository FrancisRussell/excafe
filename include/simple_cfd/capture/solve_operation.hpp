#ifndef SIMPLE_CFD_CAPTURE_SOLVE_OPERATION_HPP
#define SIMPLE_CFD_CAPTURE_SOLVE_OPERATION_HPP

#include <map>
#include "fields/named_field.hpp"
#include "fields/field.hpp"
#include "fields/discrete_expr_container.hpp"

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

    detail::DiscreteExprContainer exprContainer;
  }

  void execute()
  {
    // TODO: Implement me!
  }
};

}

#endif
