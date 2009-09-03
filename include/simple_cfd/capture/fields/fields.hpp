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

namespace cfd
{

Field linear_solve(const Operator& A, const Field& b)
{
 return Field(new detail::LinearSolve(A.getExpr(), b.getExpr()));
}

}

#endif
