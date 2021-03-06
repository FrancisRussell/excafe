#ifndef EXCAFE_CAPTURE_FIELDS_FIELDS_HPP
#define EXCAFE_CAPTURE_FIELDS_FIELDS_HPP

#include "function_space.hpp"
#include "element.hpp"
#include "scalar.hpp"
#include "field.hpp"
#include "named_field.hpp"
#include "boundary_condition.hpp"
#include "operator.hpp"
#include "temporal_index.hpp"
#include "temporal_index_expr.hpp"
#include "indexed_holder.hpp"
#include "linear_solve.hpp"
#include "temporal_index_offset.hpp"
#include "discrete_field_projection.hpp"
#include "linear_system.hpp"

namespace excafe
{

namespace detail
{

struct final_tag
{
  operator detail::TemporalIndexOffset() const
  {
    return detail::TemporalIndexOffset(detail::TemporalIndexOffset::final_tag(), 0);
  }
};

}

detail::final_tag final;

Field linear_solve(const Operator& A, const Field& b)
{
  const FunctionSpace trialSpace(A.getExpr()->getTrialSpace());
  const Field guess(new detail::DiscreteFieldZero(trialSpace));
  return Field(new detail::LinearSolve(A.getExpr(), guess.getExpr(), b.getExpr()));
}

detail::TemporalIndexExpr operator-(const TemporalIndex& e, const signed offset)
{
  return detail::TemporalIndexExpr::relative(e.getIndex(), -offset);
}

detail::TemporalIndexOffset operator-(const detail::final_tag&, const signed offset)
{
  return detail::TemporalIndexOffset(detail::TemporalIndexOffset::final_tag(), -offset);
}

Field project(const Field& field, const FunctionSpace& functionSpace)
{
  return Field(new detail::DiscreteFieldProjection(field.getExpr(), functionSpace.getExpr()));
}

LinearSystem assembleGalerkinSystem(const FunctionSpace& functionSpace, const forms::BilinearFormIntegralSum& lhs,
                                    const Field& rhs, const BoundaryCondition& bc, const Field& initialGuess)
{
  return LinearSystem(functionSpace.getExpr(), functionSpace.getExpr(), lhs, initialGuess.getExpr(), rhs.getExpr(), bc);
}

LinearSystem assembleGalerkinSystem(const FunctionSpace& functionSpace, const forms::BilinearFormIntegralSum& lhs,
                                    const Field& rhs, const BoundaryCondition& bc)
{
  const Field initialGuess(new detail::DiscreteFieldZero(functionSpace));
  return assembleGalerkinSystem(functionSpace, lhs, rhs, bc, initialGuess);
}

}

#endif
