#ifndef SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_FRACTION_HPP
#define SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_FRACTION_HPP

#include "numeric_fwd.hpp"
#include <numeric/polynomial_fraction.hpp>
#include <numeric/optimised_polynomial.hpp>
#include <vector>
#include <cstddef>
#include <cassert>

namespace cfd
{

template<typename V>
class OptimisedPolynomialFraction
{
public:
  typedef V variable_t;

private:
  typedef OptimisedPolynomial<variable_t> polynomial_t;
  polynomial_t dividend;
  polynomial_t divisor;

  double evaluate(const std::vector<double>& params) const
  {
    return dividend.evaluate(params) / divisor.evaluate(params);
  }

public:
  OptimisedPolynomialFraction(const PolynomialFraction<variable_t>& p) : dividend(p.getDividend()),
    divisor(p.getDivisor())
  {

  }
};

}

#endif
