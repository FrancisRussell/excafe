#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_APPLY_BC_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_APPLY_BC_HPP

#include <set>
#include <memory>
#include "operator_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include <simple_cfd/capture/indices/propagation_rule.hpp>
#include <simple_cfd/capture/indices/propagation_rules.hpp>
#include <simple_cfd/capture/indices/index_propagation_all.hpp>

namespace cfd
{

namespace detail
{

class OperatorApplyBC : public OperatorExpr
{
private:
  OperatorExpr::expr_ptr op;
  boost::any bc;

public:
  OperatorApplyBC(const OperatorExpr::expr_ptr& _op, const boost::any& _bc) : op(_op), bc(_bc)
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

  boost::any getBoundaryCondition() const
  {
    return bc;
  }
};


}

}

#endif
