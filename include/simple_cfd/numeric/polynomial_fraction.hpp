#ifndef SIMPLE_CFD_NUMERIC_POLYNOMIAL_FRACTION_HPP
#define SIMPLE_CFD_NUMERIC_POLYNOMIAL_FRACTION_HPP

#include <cstddef>
#include <boost/operators.hpp>
#include "polynomial.hpp"

namespace cfd
{

template<typename V>
class PolynomialFraction : boost::addable<PolynomialFraction<V>, double,
                           boost::subtractable<PolynomialFraction<V>, double,
                           boost::dividable<PolynomialFraction<V>, double,
                           boost::multipliable<PolynomialFraction<V>, double, 
                           boost::addable<PolynomialFraction<V>,
                           boost::subtractable< PolynomialFraction<V>,
                           boost::dividable< PolynomialFraction<V>,
                           boost::multipliable< PolynomialFraction<V>
                           > > > > > > > >
{
public:
  typedef V variable_t;
  typedef Polynomial<variable_t> polynomial_t;
  typedef OptimisedPolynomialFraction<variable_t> optimised_t;

private:
  polynomial_t dividend;
  polynomial_t divisor;

  PolynomialFraction(const polynomial_t& _dividend, const polynomial_t& _divisor) :
    dividend(_dividend), divisor(_divisor)
  {
  }

public:
  PolynomialFraction() : dividend(0.0), divisor(1.0)
  {
  }

  PolynomialFraction(const double constant) : dividend(constant), divisor(1.0)
  {
  }

  PolynomialFraction(const polynomial_t& p) : dividend(p), divisor(1.0)
  {
  }

  PolynomialFraction(const variable_t& variable) : dividend(variable), divisor(1.0)
  {
  }

  PolynomialFraction(const variable_t& variable, const std::size_t exponent) :
    dividend(variable, exponent), divisor(1.0)
  {
  }

  PolynomialFraction(const double coefficient, const variable_t& variable, const std::size_t exponent) :
    dividend(coefficient, variable, exponent), divisor(1.0)
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

  void replaceIndependentVariable(const variable_t& from, const variable_t& to)
  {
    dividend.replaceIndependentVariable(from, to);
    divisor.replaceIndependentVariable(from, to);
  }

  PolynomialFraction substituteValue(const variable_t& variable, const double value) const
  {
    polynomial_t newDividend = dividend.substituteValue(variable, value);
    polynomial_t newDivisor = divisor.substituteValue(variable, value);
    return PolynomialFraction(newDividend, newDivisor);
  }

  polynomial_t getDividend() const
  {
    return dividend;
  }

  polynomial_t getDivisor() const
  {
    return divisor;
  }

  PolynomialFraction& operator*=(const double x)
  {
    dividend *= x;
    return *this;
  }

  PolynomialFraction& operator/=(const double x)
  {
    divisor *= x;
    return *this;
  }

  PolynomialFraction& operator+=(const double x)
  {
    dividend += divisor*x;
    return *this;
  }

  PolynomialFraction& operator-=(const double x)
  {
    dividend -= divisor*x;
    return *this;
  }

  PolynomialFraction& operator*=(const PolynomialFraction& p)
  {
    dividend *= p.dividend;
    divisor *= p.divisor;
    return *this;
  }

  PolynomialFraction& operator/=(const PolynomialFraction& p)
  {
    dividend *= p.divisor;
    divisor *= p.dividend;
    return *this;
  }

  PolynomialFraction& operator+=(const PolynomialFraction& p)
  {
    dividend = dividend * p.divisor + p.dividend * divisor;
    divisor *= p.divisor;
    return *this;
  }

  PolynomialFraction& operator-=(const PolynomialFraction& p)
  {
    dividend = dividend * p.divisor - p.dividend * divisor;
    divisor *= p.divisor;
    return *this;
  }

  PolynomialFraction operator-() const
  {
    PolynomialFraction result(*this);
    result *= 1.0;
    return result;
  }

  optimised_t optimise() const
  {
    return optimised_t(*this);
  }

  PolynomialFraction derivative(const variable_t& variable) const
  {
    polynomial_t newDividend = dividend.derivative(variable)*divisor - dividend*divisor.derivative(variable);
    polynomial_t newDivisor = divisor * divisor;
    return PolynomialFraction(newDividend, newDivisor);
  }
};

}
#endif
