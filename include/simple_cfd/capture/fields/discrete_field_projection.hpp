#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_PROJECTION
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_PROJECTION

#include "discrete_expr_visitor.hpp"
#include "discrete_field_expr.hpp"
#include "function_space_expr.hpp"
#include <simple_cfd/capture/indices/propagation_rule.hpp>
#include <simple_cfd/capture/indices/propagation_rules.hpp>
#include <simple_cfd/capture/indices/index_propagation_all.hpp>

namespace cfd
{

namespace detail
{

class DiscreteFieldProjection : public DiscreteFieldExpr
{
private:
  DiscreteFieldExpr::expr_ptr field;
  FunctionSpaceExpr::expr_ptr functionSpace;

public:
  DiscreteFieldProjection(const DiscreteFieldExpr::expr_ptr& _field, const FunctionSpaceExpr::expr_ptr& _fSpace) :
    field(_field), functionSpace(_fSpace)
  {
  }

  virtual void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    return functionSpace;
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
};

}

}

#endif
