#include <simple_cfd/capture/fields/scalar.hpp>
#include <simple_cfd/capture/fields/scalar_expr.hpp>
#include <simple_cfd/capture/fields/scalar_literal.hpp>
#include <simple_cfd/capture/fields/scalar_binary_operator.hpp>

namespace cfd
{

Scalar::Scalar() : expr(new detail::ScalarLiteral(0))
{
}

Scalar::Scalar(const double s) : expr(new detail::ScalarLiteral(s))
{
}

Scalar::Scalar(detail::ScalarExpr* const _expr) : expr(_expr)
{
}

Scalar::Scalar(const expr_ptr _expr) : expr(_expr)
{
}

Scalar& Scalar::operator+=(const Scalar& s)
{
  expr = expr_ptr(new detail::ScalarBinaryOperator(expr, s.getExpr(), detail::ScalarBinaryOperator::add_tag()));
  return *this;
}

Scalar& Scalar::operator-=(const Scalar& s)
{
  expr = expr_ptr(new detail::ScalarBinaryOperator(expr, s.getExpr(), detail::ScalarBinaryOperator::sub_tag()));
  return *this;
}

Scalar& Scalar::operator*=(const Scalar& s)
{
  expr = expr_ptr(new detail::ScalarBinaryOperator(expr, s.getExpr(), detail::ScalarBinaryOperator::mul_tag()));
  return *this;
}

Scalar& Scalar::operator/=(const Scalar& s)
{
  expr = expr_ptr(new detail::ScalarBinaryOperator(expr, s.getExpr(), detail::ScalarBinaryOperator::div_tag()));
  return *this;
}

Scalar Scalar::operator-() const
{
  return -1.0 * (*this);
}

Scalar& Scalar::operator=(const Scalar& s)
{
  expr = s.expr;
  return *this;
}

Scalar Scalar::operator<(const Scalar& s)
{
  return Scalar(new detail::ScalarBinaryOperator(expr, s.getExpr(), detail::ScalarBinaryOperator::lt_tag()));
}

Scalar Scalar::operator<=(const Scalar& s)
{
  return Scalar(new detail::ScalarBinaryOperator(expr, s.getExpr(), detail::ScalarBinaryOperator::lte_tag()));
}

Scalar Scalar::operator>(const Scalar& s)
{
  return Scalar(new detail::ScalarBinaryOperator(expr, s.getExpr(), detail::ScalarBinaryOperator::gt_tag()));
}

Scalar Scalar::operator>=(const Scalar& s)
{
  return Scalar(new detail::ScalarBinaryOperator(expr, s.getExpr(), detail::ScalarBinaryOperator::gte_tag()));
}

Scalar Scalar::operator==(const Scalar& s)
{
  return Scalar(new detail::ScalarBinaryOperator(expr, s.getExpr(), detail::ScalarBinaryOperator::eq_tag()));
}

Scalar::expr_ptr Scalar::getExpr() const
{
  return expr;
}

}
