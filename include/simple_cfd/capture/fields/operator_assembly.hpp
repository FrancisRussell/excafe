#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ASSEMBLY_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ASSEMBLY_HPP

#include "operator_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include <simple_cfd/capture/forms/bilinear_form_integral_sum.hpp>

namespace cfd
{

namespace detail
{

class OperatorAssembly : public OperatorExpr
{
private:
  const FunctionSpace trialSpace;
  const FunctionSpace testSpace;
  const forms::BilinearFormIntegralSum sum;

public:
  OperatorAssembly(const FunctionSpace& _trialSpace, const FunctionSpace& _testSpace,
                   const forms::BilinearFormIntegralSum& _sum) : 
    trialSpace(_trialSpace), testSpace(_testSpace), sum(_sum)
  {
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
