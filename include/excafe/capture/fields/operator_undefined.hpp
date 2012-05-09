#ifndef EXCAFE_CAPTURE_FIELDS_OPERATOR_UNDEFINED_HPP
#define EXCAFE_CAPTURE_FIELDS_OPERATOR_UNDEFINED_HPP

#include "function_space_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_undefined.hpp"
#include "operator_expr.hpp"
#include "temporal_index_set.hpp"
#include <excafe/capture/indices/propagation_rules.hpp>

namespace excafe
{

namespace detail
{

class OperatorUndefined : public OperatorExpr
{
public:
  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getTrialSpace() const
  {
    return FunctionSpaceExpr::expr_ptr(new FunctionSpaceUndefined());
  }

  virtual FunctionSpaceExpr::expr_ptr getTestSpace() const
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
