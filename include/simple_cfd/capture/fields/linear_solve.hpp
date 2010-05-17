#ifndef SIMPLE_CFD_CAPTURE_FIELDS_LINEAR_SOLVE_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_LINEAR_SOLVE_HPP

#include <cassert>
#include "discrete_expr_visitor.hpp"
#include "discrete_field_expr.hpp"
#include "temporal_index_set.hpp"
#include <simple_cfd/capture/indices/propagation_rule.hpp>
#include <simple_cfd/capture/indices/propagation_rules.hpp>
#include <simple_cfd/capture/indices/index_propagation_all.hpp>

namespace cfd
{

namespace detail
{

class LinearSolve : public DiscreteFieldExpr
{
private:
  OperatorExpr::expr_ptr operation;
  DiscreteFieldExpr::expr_ptr initialGuess;
  DiscreteFieldExpr::expr_ptr operand;

public:
  LinearSolve(const OperatorExpr::expr_ptr& _operation, const DiscreteFieldExpr::expr_ptr& _initialGuess,
    const DiscreteFieldExpr::expr_ptr& _operand) :
    operation(_operation), initialGuess(_initialGuess), operand(_operand)
  {
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    assert(operation->getTestSpace() == operand->getFunctionSpace());
    return operation->getTrialSpace();
  }

  OperatorExpr& getOperator() const
  {
    return *operation;
  }

  DiscreteFieldExpr& getField() const
  {
    return *operand;
  }

  DiscreteFieldExpr& getInitialGuess() const
  {
    return *initialGuess;
  }

  virtual PropagationRules getPropagationRules()
  {
    PropagationRules rules;
    rules.insert(std::auto_ptr<PropagationRule>(new IndexPropagationAll(*operation, *this)));
    rules.insert(std::auto_ptr<PropagationRule>(new IndexPropagationAll(*initialGuess, *this)));
    rules.insert(std::auto_ptr<PropagationRule>(new IndexPropagationAll(*operand, *this)));
    return rules;
  }

  virtual std::set<DiscreteExpr*> getDependencies() const 
  {
    std::set<DiscreteExpr*> dependencies;
    dependencies.insert(&(*operation));
    dependencies.insert(&(*initialGuess));
    dependencies.insert(&(*operand));
    return dependencies;
  }
};

}

}

#endif
