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

  void replaceIndependentVariable(const variable_t& from, const variable_t& to)
  {
    dividend.replaceIndependentVariable(from, to);
    divisor.replaceIndependentVariable(from, to);
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

  PolynomialFraction derivative(const variable_t& variable) const
  {
    polynomial_t newDividend = dividend.derivative(variable)*divisor - dividend*divisor.derivative(variable);
    polynomial_t newDivisor = divisor * divisor;
    return PolynomialFraction(newDividend, newDivisor);
  }
};

}
#endif
