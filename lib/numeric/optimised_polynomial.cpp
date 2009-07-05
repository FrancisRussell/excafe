#include <numeric/optimised_polynomial.hpp>
#include <numeric/polynomial.hpp>
#include <numeric/monomial.hpp>
#include <utility>
#include <vector>
#include <set>
#include <cstddef>
#include <cassert>
#include <numeric>

namespace cfd
{

namespace 
{

struct Pow
{
  inline double operator()(const double value, const std::size_t exponent) const
  {
    double result=1.0;

    for(std::size_t power=0; power<exponent; ++power)
      result *= value;

    return result;
  }
};

}

std::vector<std::size_t> OptimisedPolynomial::buildExponentVector(const Monomial& m) const
{
  std::vector<std::size_t> exponents;

  for(std::set<std::string>::const_iterator varIter(variables.begin()); varIter!=variables.end(); ++varIter)
    exponents.push_back(m.getExponent(*varIter));

  return exponents;
}

OptimisedPolynomial::OptimisedPolynomial(const Polynomial& p) : variables(p.getIndependentVariables()),
  paramData(variables.size())
{
  for(Polynomial::const_iterator mIter(p.begin()); mIter!=p.end(); ++mIter)
    coefficients.push_back(std::make_pair(buildExponentVector(mIter->first), mIter->second));
}

double OptimisedPolynomial::evaluate(const std::vector<double>& params) const
{
  double result = 0.0;

  for(coefficient_vec_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
    result += std::inner_product(params.begin(), params.end(), cIter->first.begin(), cIter->second, std::multiplies<double>(), Pow());

  return result;
}

double OptimisedPolynomial::operator()() const
{
  assert(variables.empty());
  return evaluate(paramData);
}

double OptimisedPolynomial::operator()(const double a) const
{
  assert(variables.size() == 1);

  paramData[0] = a;

  return evaluate(paramData);
}

double OptimisedPolynomial::operator()(const double a, const double b) const
{
  assert(variables.size() == 2);

  paramData[0] = a;
  paramData[1] = b;

  return evaluate(paramData);
}

double OptimisedPolynomial::operator()(const double a, const double b, const double c) const
{
  assert(variables.size() == 3);

  paramData[0] = a;
  paramData[1] = b;
  paramData[2] = c;

  return evaluate(paramData);
}

}
