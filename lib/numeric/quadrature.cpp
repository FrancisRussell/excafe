#include <simple_cfd/numeric/quadrature.hpp>
#include <simple_cfd/numeric/math_utilities.hpp>
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/numeric/optimised_polynomial.hpp>
#include <simple_cfd/numeric/cast.hpp>
#include <map>
#include <set>
#include <string>
#include <cstddef>
#include <cmath>

namespace cfd
{

std::map<double, double> Quadrature::getGauss(const std::size_t n)
{
  const std::size_t q = cfd::numeric_cast<std::size_t>(std::ceil((n + 1)/2.0));
  const std::set<double> roots = MathUtilities::jacobi_roots(0.0, 0.0, q);

  const OptimisedPolynomial<std::string> legendrePrime = MathUtilities::jacobi(0.0, 0.0, q).derivative("x").optimise();
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
  const std::size_t q = cfd::numeric_cast<std::size_t>(std::ceil((n + 2)/2.0));
  std::set<double> roots = MathUtilities::jacobi_roots(0.0, 1.0, q-1);
  roots.insert(-1.0);

  const OptimisedPolynomial<std::string> legendre = MathUtilities::jacobi(0.0, 0.0, q-1).optimise();
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
  const std::size_t q = cfd::numeric_cast<std::size_t>(std::ceil((n + 3)/2.0));
  std::set<double> roots = MathUtilities::jacobi_roots(1.0, 1.0, q-2);
  roots.insert(-1.0);
  roots.insert(1.0);

  const OptimisedPolynomial<std::string> legendre = MathUtilities::jacobi(0.0, 0.0, q-1).optimise();
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
