#include <excafe/capture/fields/operator.hpp>
#include <excafe/capture/fields/function_space.hpp>
#include <excafe/capture/fields/function_space_expr.hpp>
#include <excafe/capture/fields/operator_assembly.hpp>
#include <excafe/capture/fields/operator_application.hpp>
#include <excafe/capture/fields/operator_addition.hpp>
#include <excafe/capture/fields/operator_undefined.hpp>
#include <excafe/capture/forms/bilinear_form_integral_sum.hpp>

namespace excafe
{

Operator::Operator(const FunctionSpace& _trialSpace, const FunctionSpace& _testSpace) : trialSpace(_trialSpace.getExpr()),
  testSpace(_testSpace.getExpr()), expr(new detail::OperatorUndefined())
{
}

Operator::Operator(const detail::FunctionSpaceExpr::expr_ptr& _trialSpace, 
                   const detail::FunctionSpaceExpr::expr_ptr& _testSpace, 
                   detail::OperatorExpr* const _expr) : 
  trialSpace(_trialSpace), testSpace(_testSpace), expr(_expr)
{
}

Operator::Operator(const detail::FunctionSpaceExpr::expr_ptr& _trialSpace, 
                   const detail::FunctionSpaceExpr::expr_ptr& _testSpace, 
                   const detail::OperatorExpr::expr_ptr _expr) : 
  trialSpace(_trialSpace), testSpace(_testSpace), expr(_expr)
{
}


Operator& Operator::operator=(const forms::BilinearFormIntegralSum& sum)
{
  expr = expr_ptr(new detail::OperatorAssembly(trialSpace, testSpace, sum));
  return *this;
}

Field Operator::operator*(const Field& field)
{
  return Field(new detail::OperatorApplication(expr, field.getExpr()));
}

Operator Operator::operator+(const Operator& o) const
{
  return Operator(trialSpace, testSpace, new detail::OperatorAddition(expr, o.getExpr()));
}

Operator Operator::operator+(const forms::BilinearFormIntegralSum& sum) const
{
  const Operator assembled(trialSpace, testSpace, new detail::OperatorAssembly(trialSpace, testSpace, sum));
  return *this + assembled;
}

Operator::expr_ptr Operator::getExpr() const
{
  return expr;
}

}
