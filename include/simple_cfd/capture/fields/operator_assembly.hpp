#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ASSEMBLY_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ASSEMBLY_HPP

#include "operator_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_expr.hpp"
#include <simple_cfd/capture/forms/bilinear_form_integral_sum.hpp>

namespace cfd
{

namespace detail
{

class OperatorAssembly : public OperatorExpr
{
private:
  const FunctionSpaceExpr::expr_ptr trialSpace;
  const FunctionSpaceExpr::expr_ptr testSpace;
  const forms::BilinearFormIntegralSum sum;

public:
  OperatorAssembly(const FunctionSpace& _trialSpace, const FunctionSpace& _testSpace,
                   const forms::BilinearFormIntegralSum& _sum) : 
    trialSpace(_trialSpace.getExpr()), testSpace(_testSpace.getExpr()), sum(_sum)
  {
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getTrialSpace() const
  {
    return trialSpace;
  }

  virtual FunctionSpaceExpr::expr_ptr getTestSpace() const
  {
    return testSpace;
  }
};

}

}

#endif
