#ifndef SIMPLE_CFD_CAPTURE_FIELDS_LINEAR_SYSTEM_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_LINEAR_SYSTEM_HPP

#include "operator.hpp"
#include "boundary_condition.hpp"
#include "field.hpp"
#include "discrete_field_expr.hpp"
#include "discrete_field_apply_bc.hpp"
#include "operator_apply_bc.hpp"

namespace cfd
{

class LinearSystem
{
private:
  typedef detail::FunctionSpaceExpr::expr_ptr function_space_ptr;

  function_space_ptr trialSpace;
  function_space_ptr testSpace;
  forms::BilinearFormIntegralSum lhs;
  detail::DiscreteFieldExpr::expr_ptr rhs;
  BoundaryCondition bc;

  detail::OperatorExpr::expr_ptr constrainedSystem;
  detail::DiscreteFieldExpr::expr_ptr constrainedLoad;
  detail::DiscreteFieldExpr::expr_ptr unknown;

public:
  LinearSystem(const function_space_ptr& _trialSpace, const function_space_ptr& _testSpace, 
               const forms::BilinearFormIntegralSum& _lhs, const detail::DiscreteFieldExpr::expr_ptr& _rhs,
               const BoundaryCondition& _bc) :
    trialSpace(_trialSpace), testSpace(_testSpace), lhs(_lhs), rhs(_rhs), bc(_bc)
  {
    constrainedSystem = detail::OperatorExpr::expr_ptr(new detail::OperatorAssembly(trialSpace, testSpace, lhs));
    constrainedSystem = detail::OperatorExpr::expr_ptr(new detail::OperatorApplyBC(constrainedSystem, bc));

    constrainedLoad = rhs;
    constrainedLoad = detail::DiscreteFieldExpr::expr_ptr(new detail::DiscreteFieldApplyBC(constrainedLoad, bc));

    unknown = detail::DiscreteFieldExpr::expr_ptr(new detail::LinearSolve(constrainedSystem, constrainedLoad));
  }

  Operator getConstrainedSystem() const
  {
    return Operator(constrainedSystem->getTrialSpace(), constrainedSystem->getTestSpace(), constrainedSystem);
  }

  Field getConstrainedLoad() const
  {
    return Field(constrainedLoad);
  }

  Field getSolution() const
  {
    return Field(unknown);
  }
};

}

#endif
