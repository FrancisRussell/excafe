#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELDS_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELDS_HPP

#include "function_space.hpp"
#include "element.hpp"
#include "scalar.hpp"
#include "field.hpp"
#include "named_field.hpp"
#include "operator.hpp"
#include "temporal_index.hpp"
#include "indexed_holder.hpp"
#include "linear_solve.hpp"
#include "temporal_index_offset.hpp"

namespace cfd
{

namespace detail
{

struct final_tag {};

}

detail::final_tag final;

Field linear_solve(const Operator& A, const Field& b)
{
 return Field(new detail::LinearSolve(A.getExpr(), b.getExpr()));
}

detail::TemporalIndexExpr operator-(const TemporalIndex& e, const unsigned offset)
{
  return detail::TemporalIndexExpr::relative(e.getIndex(), offset);
}

detail::TemporalIndexOffset operator-(const detail::final_tag&, const unsigned offset)
{
  return detail::TemporalIndexOffset(detail::TemporalIndexOffset::final_tag(), offset);
}


}

#endif
