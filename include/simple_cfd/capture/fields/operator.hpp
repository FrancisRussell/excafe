#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_HPP

#include "fields_fwd.hpp"
#include "operator_expr.hpp"
#include "function_space_expr.hpp"
#include <simple_cfd/capture/forms/bilinear_form_integral_sum.hpp>

namespace cfd
{

class Operator
{
public:
  typedef detail::OperatorExpr::expr_ptr expr_ptr;

private:
  detail::FunctionSpaceExpr::expr_ptr trialSpace;
  detail::FunctionSpaceExpr::expr_ptr testSpace;
  expr_ptr expr;

public:
  Operator(const FunctionSpace& _trialSpace, const FunctionSpace& _testSpace);
  Operator(const detail::FunctionSpaceExpr::expr_ptr& _trialSpace, 
           const detail::FunctionSpaceExpr::expr_ptr& _testSpace, 
           detail::OperatorExpr* const _expr);
  Operator(const detail::FunctionSpaceExpr::expr_ptr& _trialSpace, 
           const detail::FunctionSpaceExpr::expr_ptr& _testSpace, 
           const detail::OperatorExpr::expr_ptr _expr);
  Operator& operator=(const forms::BilinearFormIntegralSum& sum);
  Operator operator+(const Operator& o) const;
  Operator operator+(const forms::BilinearFormIntegralSum& sum) const;
  Field operator*(const Field& field);
  expr_ptr getExpr() const;
};

}

#endif
