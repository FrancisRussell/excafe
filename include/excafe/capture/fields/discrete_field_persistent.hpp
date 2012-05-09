#ifndef EXCAFE_CAPTURE_FIELDS_DISCRETE_FIELD_PERSISTENT_HPP
#define EXCAFE_CAPTURE_FIELDS_DISCRETE_FIELD_PERSISTENT_HPP

#include <set>
#include <string>
#include "discrete_field_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_expr.hpp"
#include "temporal_index_set.hpp"
#include <excafe/capture/indices/propagation_rules.hpp>

namespace excafe
{

namespace detail
{

class DiscreteFieldPersistent : public DiscreteFieldExpr
{
private:
  const std::string name;
  const FunctionSpaceExpr::expr_ptr functionSpace;

public:
  DiscreteFieldPersistent(const std::string& _name, const FunctionSpaceExpr::expr_ptr& _functionSpace) :
    name(_name), functionSpace(_functionSpace)
  {
  }

  void accept(DiscreteExprVisitor& visitor)
  {
    visitor.visit(*this);
  }

  FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    return functionSpace;
  }

  PropagationRules getPropagationRules()
  {
    return PropagationRules();
  }

  virtual std::set<DiscreteExpr*> getDependencies() const 
  {
    return std::set<DiscreteExpr*>();
  }

  std::string getName() const
  {
    return name;
  }
};

}

}

#endif
