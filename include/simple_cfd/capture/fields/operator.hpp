#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_HPP

#include "fields_fwd.hpp"
#include "operator_expr.hpp"
#include <simple_cfd/capture/forms/bilinear_form_integral_sum.hpp>

namespace cfd
{

class Operator
{
public:
  typedef detail::OperatorExpr::expr_ptr expr_ptr;

private:
  FunctionSpace trialSpace;
  FunctionSpace testSpace;
  expr_ptr expr;

public:
  Operator(const FunctionSpace& _trialSpace, const FunctionSpace& _testSpace);
  Operator& operator=(const forms::BilinearFormIntegralSum& sum);
  Field operator*(const Field& field);
  expr_ptr getExpr() const;
};

}

#endif
