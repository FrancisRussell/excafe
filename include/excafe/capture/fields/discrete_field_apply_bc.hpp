#ifndef EXCAFE_CAPTURE_FIELDS_DISCRETE_FIELD_APPLY_BC_HPP
#define EXCAFE_CAPTURE_FIELDS_DISCRETE_FIELD_APPLY_BC_HPP

#include <set>
#include <memory>
#include "discrete_field_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_expr.hpp"
#include "boundary_condition.hpp"
#include <excafe/capture/indices/propagation_rules.hpp>
#include <excafe/capture/indices/propagation_rule.hpp>
#include <excafe/capture/indices/index_propagation_all.hpp>

namespace excafe
{

namespace detail
{

class DiscreteFieldApplyBC : public DiscreteFieldExpr
{
private:
  DiscreteFieldExpr::expr_ptr field;
  BoundaryCondition bc;

public:
  DiscreteFieldApplyBC(const DiscreteFieldExpr::expr_ptr& f, const BoundaryCondition& _bc) : field(f), bc(_bc)
  {
  }

  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    return field->getFunctionSpace();
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  DiscreteFieldExpr& getField() const
  {
    return *field;
  }

  virtual PropagationRules getPropagationRules()
  {
    PropagationRules rules;
    rules.insert(std::auto_ptr<PropagationRule>(new IndexPropagationAll(*field, *this)));
    return rules;
  }

  virtual std::set<DiscreteExpr*> getDependencies() const 
  {
    std::set<DiscreteExpr*> dependencies;
    dependencies.insert(&(*field));
    return dependencies;
  }

  BoundaryCondition getBoundaryCondition() const
  {
    return bc;
  }
};


}

}

#endif
