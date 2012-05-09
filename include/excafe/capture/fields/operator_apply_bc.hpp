#ifndef EXCAFE_CAPTURE_FIELDS_OPERATOR_APPLY_BC_HPP
#define EXCAFE_CAPTURE_FIELDS_OPERATOR_APPLY_BC_HPP

#include <set>
#include <memory>
#include "operator_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "boundary_condition.hpp"
#include <excafe/capture/indices/propagation_rule.hpp>
#include <excafe/capture/indices/propagation_rules.hpp>
#include <excafe/capture/indices/index_propagation_all.hpp>

namespace excafe
{

namespace detail
{

class OperatorApplyBC : public OperatorExpr
{
private:
  OperatorExpr::expr_ptr op;
  BoundaryCondition bc;

public:
  OperatorApplyBC(const OperatorExpr::expr_ptr& _op, const BoundaryCondition& _bc) : op(_op), bc(_bc)
  {
  }

  virtual void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getTrialSpace() const
  {
    return op->getTrialSpace();
  }

  virtual FunctionSpaceExpr::expr_ptr getTestSpace() const
  {
    return op->getTestSpace();
  }

  OperatorExpr& getOperator() const
  {
    return *op;
  }

  virtual PropagationRules getPropagationRules()
  {
    PropagationRules rules;
    rules.insert(std::auto_ptr<PropagationRule>(new IndexPropagationAll(*op, *this)));
    return rules;
  }

  virtual std::set<DiscreteExpr*> getDependencies() const 
  {
    std::set<DiscreteExpr*> dependencies;
    dependencies.insert(&(*op));
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
