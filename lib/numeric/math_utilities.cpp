#include <simple_cfd/numeric/math_utilities.hpp>
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/numeric/optimised_polynomial.hpp>
#include <cmath>
#include <set>
#include <vector>

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

std::set<double> MathUtilities::jacobi_roots(const double alpha, const double beta, const std::size_t n, const double epsilon)
{
  const OptimisedPolynomial j = jacobi(alpha, beta, n).optimise();
  const OptimisedPolynomial jPrime = jacobi(alpha, beta, n).derivative("x").optimise();

  std::vector<double> roots(n);

  for(std::size_t k=0; k<n; ++k)
  {
    double r = - std::cos(M_PI * (2.0*k + 1.0)/(2.0*n));

    if (k>0)
      r = (r + roots[k-1])/2.0;

    while(true)
    {
      double s = 0.0;

      for(std::size_t i=0; i<k; ++i)
        s += 1.0/(r - roots[i]);

      const double delta = -j(r) / (jPrime(r) - j(r)*s);
      r += delta;

      if (delta < epsilon)
        break;
    }

    roots[k] = r;
  }

  return std::set<double>(roots.begin(), roots.end());
}

}
