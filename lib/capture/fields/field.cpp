#include <simple_cfd/capture/fields/field.hpp>
#include <simple_cfd/capture/fields/discrete_field_undefined.hpp>
#include <simple_cfd/capture/fields/discrete_field_zero.hpp>
#include <simple_cfd/capture/fields/discrete_field_expr.hpp>
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

Field::expr_ptr Field::getExpr() const
{
  return expr;
}

}
