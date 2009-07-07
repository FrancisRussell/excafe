#ifndef SIMPLE_CFD_NUMERIC_MATH_UTILITIES_HPP
#define SIMPLE_CFD_NUMERIC_MATH_UTILITIES_HPP

#include <set>
#include <cstddef>
#include <simple_cfd/simple_cfd_fwd.hpp>

namespace cfd
{

class MathUtilities
{
private:
  static double jacobi_a_1_n(const double alpha, const double beta, const std::size_t n);
  static double jacobi_a_2_n(const double alpha, const double beta, const std::size_t n);
  static double jacobi_a_3_n(const double alpha, const double beta, const std::size_t n);
  static double jacobi_a_4_n(const double alpha, const double beta, const std::size_t n);

public:
  static Polynomial jacobi(const double alpha, const double beta, const std::size_t n);
  static std::set<double> jacobi_roots(const double alpha, const double beta, const std::size_t n, const double epsilon = 1e-8);
};

}

#endif
