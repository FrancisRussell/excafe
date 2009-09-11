#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_ELEMENT_WISE_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_ELEMENT_WISE_HPP

#include "discrete_field_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_expr.hpp"
#include <boost/variant.hpp>

namespace cfd
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

  operator_t getOperator() const
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
};

}

}
#endif
