#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_ASSEMBLY_POLYNOMIAL_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_ASSEMBLY_POLYNOMIAL_HPP

#include "scalar_placeholder.hpp"
#include <simple_cfd/numeric/polynomial_fraction.hpp>
#include <simple_cfd/numeric/optimised_polynomial_fraction.hpp>

namespace cfd
{

namespace detail
{

// typedef for easy alteration of capture polynomial type
typedef PolynomialFraction<ScalarPlaceholder> assembly_polynomial_t;
typedef OptimisedPolynomialFraction<ScalarPlaceholder> optimised_assembly_polynomial_t;

}

}

#endif
