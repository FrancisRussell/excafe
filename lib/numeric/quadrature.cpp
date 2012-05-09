#include <excafe/numeric/quadrature.hpp>
#include <excafe/numeric/math_utilities.hpp>
#include <excafe/numeric/excafe_expression.hpp>
#include <excafe/numeric/cast.hpp>
#include <map>
#include <set>
#include <string>
#include <cstddef>
#include <cmath>

namespace excafe
{

std::map<double, double> Quadrature::getGauss(const std::size_t n)
{
  const std::string x("x");
  const std::size_t q = excafe::numeric_cast<std::size_t>(std::ceil((n + 1)/2.0));
  const std::set<double> roots = MathUtilities::jacobi_roots(0.0, 0.0, q);

  const ExcafeExpression<std::string>::optimised_t legendrePrime = MathUtilities::jacobi(x, 0.0, 0.0, q).derivative(x).optimise();
  std::map<double, double> quadrature;

  for(std::set<double>::const_iterator rootIter(roots.begin()); rootIter!=roots.end(); ++rootIter)
  {
    const double xi = *rootIter;
    const double legendrePrimeVal = legendrePrime(xi);
    quadrature[xi] = 2.0 / ((1.0 - xi*xi) * legendrePrimeVal * legendrePrimeVal);
  }

  return quadrature;
}

std::map<double, double> Quadrature::getGaussRadau(const std::size_t n)
{
  const std::size_t q = excafe::numeric_cast<std::size_t>(std::ceil((n + 2)/2.0));
  std::set<double> roots = MathUtilities::jacobi_roots(0.0, 1.0, q-1);
  roots.insert(-1.0);

  const std::string x("x");
  const ExcafeExpression<std::string>::optimised_t legendre = MathUtilities::jacobi(x, 0.0, 0.0, q-1).optimise();
  std::map<double, double> quadrature;

  for(std::set<double>::const_iterator rootIter(roots.begin()); rootIter!=roots.end(); ++rootIter)
  {
    const double xi = *rootIter;
    const double legendreVal = legendre(xi);
    quadrature[xi] = (1.0 - xi) / (q*q*legendreVal*legendreVal);
  }

  return quadrature;
}

std::map<double, double> Quadrature::getGaussLobatto(const std::size_t n)
{
  const std::size_t q = excafe::numeric_cast<std::size_t>(std::ceil((n + 3)/2.0));
  std::set<double> roots = MathUtilities::jacobi_roots(1.0, 1.0, q-2);
  roots.insert(-1.0);
  roots.insert(1.0);

  const std::string x("x");
  const ExcafeExpression<std::string>::optimised_t legendre = MathUtilities::jacobi(x, 0.0, 0.0, q-1).optimise();
  std::map<double, double> quadrature;

  for(std::set<double>::const_iterator rootIter(roots.begin()); rootIter!=roots.end(); ++rootIter)
  {
    const double xi = *rootIter;
    const double legendreVal = legendre(xi);
    quadrature[xi] = 2.0 / (q*(q-1)*legendreVal*legendreVal);
  }

  return quadrature;
}

}
