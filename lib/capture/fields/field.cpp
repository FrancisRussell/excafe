#include <simple_cfd/capture/fields/field.hpp>
#include <simple_cfd/capture/fields/discrete_field_undefined.hpp>
#include <simple_cfd/capture/fields/discrete_field_zero.hpp>
#include <simple_cfd/capture/fields/discrete_field_expr.hpp>
#include <simple_cfd/capture/fields/discrete_field_element_wise.hpp>
#include <simple_cfd/capture/fields/discrete_field_two_norm.hpp>
#include <simple_cfd/capture/fields/operator.hpp>

namespace cfd
{

Field::Field() : expr(new detail::DiscreteFieldUndefined())
{
}

Field::Field(const FunctionSpace& functionSpace) : expr(new detail::DiscreteFieldZero(functionSpace))
{
}

Field::Field(detail::DiscreteFieldExpr* const _expr) : expr(_expr)
{
}

Field::Field(const expr_ptr _expr) : expr(_expr)
{
}

Field& Field::operator+=(const Field& f)
{
  expr = expr_ptr(new detail::DiscreteFieldElementWise(expr, f.expr, detail::DiscreteFieldElementWise::add_tag()));
  return *this;
}

Field& Field::operator-=(const Field& f)
{
  expr = expr_ptr(new detail::DiscreteFieldElementWise(expr, f.expr, detail::DiscreteFieldElementWise::sub_tag()));
  return *this;
}

Scalar Field::two_norm() const
{
  return Scalar(new detail::DiscreteFieldTwoNorm(expr));
}

Field::expr_ptr Field::getExpr() const
{
  return expr;
}

}
