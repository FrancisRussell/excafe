#ifndef EXCAFE_CAPTURE_FIELDS_SCALAR_BINARY_OPERATOR_HPP
#define EXCAFE_CAPTURE_FIELDS_SCALAR_BINARY_OPERATOR_HPP

#include "scalar_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "temporal_index_set.hpp"
#include <boost/variant.hpp>
#include <excafe/capture/indices/propagation_rules.hpp>
#include <excafe/capture/indices/propagation_rule.hpp>
#include <excafe/capture/indices/index_propagation_all.hpp>

namespace excafe
{

namespace detail
{

class ScalarBinaryOperator : public ScalarExpr
{
public:
  class add_tag {};
  class sub_tag {};
  class div_tag {};
  class mul_tag {};
  class lt_tag {};
  class gt_tag {};
  class lte_tag {};
  class gte_tag {};
  class eq_tag {};
  typedef boost::variant<add_tag,sub_tag,div_tag,mul_tag,lt_tag,gt_tag,lte_tag,gte_tag,eq_tag> operator_t;

private:
  operator_t operation;
  ScalarExpr::expr_ptr left;
  ScalarExpr::expr_ptr right;

public:
  template<typename operator_type>
  ScalarBinaryOperator(const ScalarExpr::expr_ptr& _left, const ScalarExpr::expr_ptr& _right, const operator_type operator_tag) : 
    operation(operator_tag), left(_left), right(_right)
  {
  }

  ScalarExpr& getLeft() const
  {
    return *left;
  }

  ScalarExpr& getRight() const
  {
    return *right;
  }

  operator_t getOperator() const
  {
    return operation;
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
