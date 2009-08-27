#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_APPLICATION_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_APPLICATION_HPP

#include "discrete_field_expr.hpp"
#include "discrete_expr_visitor.hpp"

namespace cfd
{

namespace detail
{

class OperatorApplication : public DiscreteFieldExpr
{
public:
  typedef OperatorExpr::expr_ptr operator_ptr;
  typedef DiscreteFieldExpr::expr_ptr field_ptr;

private:
  operator_ptr operation;
  field_ptr field;

public:
  OperatorApplication(const operator_ptr& o, const field_ptr& f) : operation(o), field(f)
  {
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.enter(*this);
    operation->accept(v);
    field->accept(v);
    v.exit(*this);
  }
};

}

}

#endif
