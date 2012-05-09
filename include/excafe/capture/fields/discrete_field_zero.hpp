#ifndef EXCAFE_CAPTURE_FIELDS_DISCRETE_FIELD_ZERO_HPP
#define EXCAFE_CAPTURE_FIELDS_DISCRETE_FIELD_ZERO_HPP

#include "function_space.hpp"
#include "function_space_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "discrete_field_expr.hpp"
#include "temporal_index_set.hpp"
#include <excafe/capture/indices/propagation_rules.hpp>

namespace excafe
{

namespace detail
{

class DiscreteFieldZero : public DiscreteFieldExpr
{
private:
  const FunctionSpaceExpr::expr_ptr functionSpace;

public:
  DiscreteFieldZero(const FunctionSpace& _functionSpace) : functionSpace(_functionSpace.getExpr())
  {
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    return functionSpace;
  }

  virtual PropagationRules getPropagationRules()
  {
    return PropagationRules();
  }

  virtual std::set<DiscreteExpr*> getDependencies() const 
  {
    return std::set<DiscreteExpr*>();
  }
};

}

}

#endif
