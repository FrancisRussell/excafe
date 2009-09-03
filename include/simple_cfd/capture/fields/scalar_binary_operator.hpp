#ifndef SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_BINARY_OPERATOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_BINARY_OPERATOR_HPP

#include "scalar_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include <boost/variant.hpp>

namespace cfd
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
    v.enter(*this);
    left->accept(v);
    right->accept(v);
    v.exit(*this);
  }
};

}

}

#endif
