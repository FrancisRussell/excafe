#ifndef SIMPLE_CFD_CAPTURE_DIMENSIONLESS_SCENARIO_HPP
#define SIMPLE_CFD_CAPTURE_DIMENSIONLESS_SCENARIO_HPP

#include "capture_fwd.hpp"
#include "fields/function_space_expr.hpp"

namespace cfd
{

namespace detail
{

class DimensionlessScenario
{
public:
  virtual void resolveFunctionSpaces(const std::set<FunctionSpaceExpr*> functionSpaces) = 0;
  virtual void optimiseLocalAssemblies(const DiscreteExprSet<discrete_operator_tag>& operators) = 0;
  virtual void execute(SolveOperation& o) = 0;
  virtual ~DimensionlessScenario() {}
};

}

}

#endif
