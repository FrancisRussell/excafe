#ifndef SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_FRACTION_HPP
#define SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_FRACTION_HPP

#include <set>
#include <vector>
#include <map>
#include <cstddef>
#include <cassert>
#include "numeric_fwd.hpp"
#include "polynomial_fraction.hpp"
#include "optimised_polynomial.hpp"

namespace cfd
{

template<typename V>
class OptimisedPolynomialFraction
{
public:
  typedef V variable_t;
  typedef double value_type;

private:
  typedef OptimisedPolynomial<variable_t> polynomial_t;
  polynomial_t dividend;
  polynomial_t divisor;

public:
  typedef typename polynomial_t::value_map value_map;

  OptimisedPolynomialFraction()
  {
  }

  OptimisedPolynomialFraction(const PolynomialFraction<variable_t>& p) : dividend(p.getDividend()),
    divisor(p.getDivisor())
  {

  }

  std::set<variable_t> getVariables() const
  {
    const std::set<variable_t> dividendVariables(dividend.getVariables());
    const std::set<variable_t> divisorVariables(divisor.getVariables());

    std::set<variable_t> result;
    result.insert(dividendVariables.begin(), dividendVariables.end());
    result.insert(divisorVariables.begin(), divisorVariables.end());
    return result;
  }

  value_type evaluate(const value_map& variableValues) const
  {
    return dividend.evaluate(variableValues)/divisor.evaluate(variableValues);
  }
};

}

#endif
