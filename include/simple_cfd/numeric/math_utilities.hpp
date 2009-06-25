#ifndef SIMPLE_CFD_NUMERIC_MATH_UTILITIES_HPP
#define SIMPLE_CFD_NUMERIC_MATH_UTILITIES_HPP

#include <cstddef>

namespace cfd
{

class MathUtilities
{
public:
  static double rising_factorial(const double x, const std::size_t n);
};

}

#endif
