#include <simple_cfd/numeric/math_utilities.hpp>

namespace cfd
{

double MathUtilities::rising_factorial(const double x, const std::size_t n)
{
  double r = 1.0;

  for(std::size_t m=0; m<n; ++m)
    r *= x+m;

  return r;
}

}
