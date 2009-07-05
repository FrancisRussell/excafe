#include <simple_cfd/numeric/math_utilities.hpp>
#include <simple_cfd/numeric/polynomial.hpp>

namespace cfd
{

double MathUtilities::rising_factorial(const double x, const std::size_t n)
{
  double r = 1.0;

  for(std::size_t m=0; m<n; ++m)
    r *= x+m;

  return r;
}

double MathUtilities::jacobi_a_1_n(const double alpha, const double beta, const std::size_t n)
{
  return 2*(n+1)*(n + alpha + beta + 1)*(2*n + alpha + beta);
}

double MathUtilities::jacobi_a_2_n(const double alpha, const double beta, const std::size_t n)
{
  return (2*n + alpha + beta + 1)*(alpha*alpha + beta*beta);
}

double MathUtilities::jacobi_a_3_n(const double alpha, const double beta, const std::size_t n)
{
  return (2*n + alpha + beta)*(2*n + alpha + beta + 1)*(2*n + alpha + beta + 2);
}

double MathUtilities::jacobi_a_4_n(const double alpha, const double beta, const std::size_t n)
{
  return 2*(n + alpha)*(n + beta)*(2*n + alpha + beta + 2);
}

Polynomial MathUtilities::jacobi(const double alpha, const double beta, const std::size_t n)
{
  if (n == 0)
  {
    return Polynomial(1.0);
  }
  else if (n == 1)
  {
    return (alpha - beta + (alpha + beta + 2.0)*Polynomial("x")) * 0.5;
  }
  else
  {
    return ((jacobi_a_2_n(alpha, beta, n-1) + jacobi_a_3_n(alpha, beta, n-1)*Polynomial("x")) * jacobi(alpha, beta, n-1) -
      jacobi_a_4_n(alpha, beta, n-1) * jacobi(alpha, beta, n-2)) / jacobi_a_1_n(alpha, beta, n-1);
  }
}

}
