#ifndef EXCAFE_CAPTURE_FIELDS_DISCRETE_FIELD_EMPTY_HPP
#define EXCAFE_CAPTURE_FIELDS_DISCRETE_FIELD_EMPTY_HPP

#include "discrete_field_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_expr.hpp"
#include "function_space_undefined.hpp"
#include "temporal_index_set.hpp"
#include <excafe/capture/indices/propagation_rules.hpp>

namespace excafe
{

namespace detail
{

class DiscreteFieldUndefined : public DiscreteFieldExpr
{
public:
  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    return FunctionSpaceExpr::expr_ptr(new FunctionSpaceUndefined());
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
