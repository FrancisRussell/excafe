#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_ASSEMBLY_POLYNOMIAL_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_ASSEMBLY_POLYNOMIAL_HPP

#include "assembly_fwd.hpp"
#include "scalar_placeholder.hpp"
#include <simple_cfd/numeric/ginac_expression.hpp>

namespace cfd
{

namespace detail
{

// typedef for easy alteration of capture polynomial type
typedef GinacExpression<ScalarPlaceholder> assembly_polynomial_t;
typedef GinacExpression<ScalarPlaceholder> optimised_assembly_polynomial_t;

}

}

#endif
