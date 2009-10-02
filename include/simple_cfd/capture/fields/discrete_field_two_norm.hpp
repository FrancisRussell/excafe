#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_TWO_NORM_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_TWO_NORM_HPP

#include "scalar_expr.hpp"
#include "discrete_field_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "temporal_index_set.hpp"
#include <simple_cfd/capture/indices/propagation_rules.hpp>
#include <simple_cfd/capture/indices/propagation_rule.hpp>
#include <simple_cfd/capture/indices/index_propagation_all.hpp>

namespace cfd
{

namespace detail
{

class DiscreteFieldTwoNorm : public ScalarExpr
{
private:
  DiscreteFieldExpr::expr_ptr field;

public:
  DiscreteFieldTwoNorm(const DiscreteFieldExpr::expr_ptr& f) : field(f)
  {
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
};

}

}

#endif
