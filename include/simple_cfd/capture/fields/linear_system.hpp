#ifndef SIMPLE_CFD_CAPTURE_FIELDS_LINEAR_SYSTEM_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_LINEAR_SYSTEM_HPP

#include <boost/function.hpp>
#include "discrete_field_expr.hpp"
#include "discrete_field_apply_bc.hpp"
#include "operator_apply_bc.hpp"

namespace cfd
{

class LinearSystem
{
private:
  typedef boost::function<void (DiscreteOperator<2>&)> system_matrix_bc_t;
  typedef boost::function<void (DiscreteField<2>&)> load_vector_bc_t;
  typedef detail::FunctionSpaceExpr::expr_ptr function_space_ptr;

  function_space_ptr trialSpace;
  function_space_ptr testSpace;
  forms::BilinearFormIntegralSum lhs;
  detail::DiscreteFieldExpr::expr_ptr rhs;
  system_matrix_bc_t systemMatrixBC;
  load_vector_bc_t loadVectorBC;

  detail::OperatorExpr::expr_ptr constrainedSystem;
  detail::DiscreteFieldExpr::expr_ptr constrainedLoad;
  detail::DiscreteFieldExpr::expr_ptr unknown;

public:
  template<typename bc_t>
  LinearSystem(const function_space_ptr& _trialSpace, const function_space_ptr& _testSpace, 
               const forms::BilinearFormIntegralSum& _lhs, const detail::DiscreteFieldExpr::expr_ptr& _rhs, const bc_t& bc) :
    trialSpace(_trialSpace), testSpace(_testSpace), lhs(_lhs), rhs(_rhs), systemMatrixBC(bc), loadVectorBC(bc)
  {
    constrainedSystem = detail::OperatorExpr::expr_ptr(new detail::OperatorAssembly(trialSpace, testSpace, lhs));
    constrainedSystem = detail::OperatorExpr::expr_ptr(new detail::OperatorApplyBC(constrainedSystem));

    constrainedLoad = rhs;
    constrainedLoad = detail::DiscreteFieldExpr::expr_ptr(new detail::DiscreteFieldApplyBC(constrainedLoad));

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
