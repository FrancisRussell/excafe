#include <simple_cfd/capture/fields/operator.hpp>
#include <simple_cfd/capture/fields/function_space.hpp>
#include <simple_cfd/capture/fields/operator_assembly.hpp>
#include <simple_cfd/capture/fields/operator_application.hpp>
#include <simple_cfd/capture/fields/operator_undefined.hpp>
#include <simple_cfd/capture/forms/bilinear_form_integral_sum.hpp>

namespace cfd
{

Operator::Operator(const FunctionSpace& _trialSpace, const FunctionSpace& _testSpace) : trialSpace(_trialSpace),
    testSpace(_testSpace), expr(new detail::OperatorUndefined())
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

}
