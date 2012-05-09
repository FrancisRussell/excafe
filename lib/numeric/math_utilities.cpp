#include <excafe/numeric/math_utilities.hpp>
#include <excafe/numeric/excafe_expression.hpp>
#include <cmath>
#include <set>
#include <vector>
#include <string>

namespace excafe
{

double MathUtilities::jacobi_a_1_n(const double alpha, const double beta, const std::size_t n)
{
  return 2*(n+1)*(n + alpha + beta + 1)*(2*n + alpha + beta);
}

double MathUtilities::jacobi_a_2_n(const double alpha, const double beta, const std::size_t n)
{
  return (2*n + alpha + beta + 1)*(alpha*alpha - beta*beta);
}

double MathUtilities::jacobi_a_3_n(const double alpha, const double beta, const std::size_t n)
{
  return (2*n + alpha + beta)*(2*n + alpha + beta + 1)*(2*n + alpha + beta + 2);
}

double MathUtilities::jacobi_a_4_n(const double alpha, const double beta, const std::size_t n)
{
  return 2*(n + alpha)*(n + beta)*(2*n + alpha + beta + 2);
}

std::set<double> MathUtilities::jacobi_roots(const double alpha, const double beta, const std::size_t n, const double epsilon)
{
  const std::string x("x");
  const ExcafeExpression<std::string>::optimised_t j = jacobi(x, alpha, beta, n).optimise();
  const ExcafeExpression<std::string>::optimised_t jPrime = jacobi(x, alpha, beta, n).derivative(x).optimise();

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
