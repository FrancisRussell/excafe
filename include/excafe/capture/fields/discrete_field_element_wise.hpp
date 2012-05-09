#ifndef EXCAFE_CAPTURE_FIELDS_DISCRETE_FIELD_ELEMENT_WISE_HPP
#define EXCAFE_CAPTURE_FIELDS_DISCRETE_FIELD_ELEMENT_WISE_HPP

#include <set>
#include <memory>
#include "discrete_field_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_expr.hpp"
#include <boost/variant.hpp>
#include <excafe/capture/indices/propagation_rules.hpp>
#include <excafe/capture/indices/propagation_rule.hpp>
#include <excafe/capture/indices/index_propagation_all.hpp>

namespace excafe
{

namespace detail
{

class DiscreteFieldElementWise : public DiscreteFieldExpr
{
public:
  class add_tag {};
  class sub_tag {};
  typedef boost::variant<add_tag,sub_tag> operator_t;

private:
  operator_t operation;
  DiscreteFieldExpr::expr_ptr left;
  DiscreteFieldExpr::expr_ptr right;

public:
  template<typename operator_type>
  DiscreteFieldElementWise(const DiscreteFieldExpr::expr_ptr& _left, const DiscreteFieldExpr::expr_ptr& _right, 
    const operator_type operator_tag) : 
    operation(operator_tag), left(_left), right(_right)
  {
  }

  DiscreteFieldExpr& getLeft() const
  {
    return *left;
  }

  DiscreteFieldExpr& getRight() const
  {
    return *right;
  }

  operator_t getOperation() const
  {
    return operation;
  }

  virtual FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    assert(left->getFunctionSpace() == right->getFunctionSpace());
    return left->getFunctionSpace();
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual PropagationRules getPropagationRules()
  {
    PropagationRules rules;
    rules.insert(std::auto_ptr<PropagationRule>(new IndexPropagationAll(*left, *this)));
    rules.insert(std::auto_ptr<PropagationRule>(new IndexPropagationAll(*right, *this)));
    return rules;
  }

  virtual std::set<DiscreteExpr*> getDependencies() const 
  {
    std::set<DiscreteExpr*> dependencies;
    dependencies.insert(&(*left));
    dependencies.insert(&(*right));
    return dependencies;
  }
};

}

}
#endif
