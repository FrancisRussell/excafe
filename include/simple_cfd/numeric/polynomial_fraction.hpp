#ifndef SIMPLE_CFD_NUMERIC_POLYNOMIAL_FRACTION_HPP
#define SIMPLE_CFD_NUMERIC_POLYNOMIAL_FRACTION_HPP

#include <cstddef>
#include <map>
#include <boost/operators.hpp>
#include "polynomial.hpp"
#include <ostream>

namespace cfd
{

template<typename V>
class PolynomialFraction : boost::arithmetic<PolynomialFraction<V>, typename Polynomial<V>::value_type,
                           boost::arithmetic<PolynomialFraction<V>,
                           boost::totally_ordered< PolynomialFraction<V>
                           > > >
{
public:
  typedef V variable_t;
  typedef Polynomial<variable_t> polynomial_t;
  typedef typename polynomial_t::value_type value_type;
  typedef typename polynomial_t::value_map value_map;
  typedef OptimisedPolynomialFraction<variable_t> optimised_t;

private:
  polynomial_t dividend;
  polynomial_t divisor;

  PolynomialFraction(const polynomial_t& _dividend, const polynomial_t& _divisor) :
    dividend(_dividend), divisor(_divisor)
  {
    simplify();
  }

  void simplify()
  {
    if (dividend == polynomial_t(0.0))
    {
      divisor = polynomial_t(1.0);
    }

    if (dividend == divisor)
    {
      dividend = polynomial_t(1.0);
      divisor = polynomial_t(1.0);
    }
  }

public:
  PolynomialFraction() : dividend(0.0), divisor(1.0)
  {
  }

  PolynomialFraction(const value_type constant) : dividend(constant), divisor(1.0)
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

  PolynomialFraction(const value_type coefficient, const variable_t& variable, const std::size_t exponent) :
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

  PolynomialFraction substituteValues(const value_map& valueMap) const
  {
    const polynomial_t newDividend = dividend.substituteValues(valueMap);
    const polynomial_t newDivisor = divisor.substituteValues(valueMap);
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

  bool operator==(const PolynomialFraction& p) const
  {
    return dividend == p.dividend && divisor == p.divisor;
  }

  PolynomialFraction& operator*=(const value_type x)
  {
    dividend *= x;
    return *this;
  }

  PolynomialFraction& operator/=(const value_type x)
  {
    dividend /= x;
    return *this;
  }

  PolynomialFraction& operator+=(const value_type x)
  {
    dividend += divisor*x;
    return *this;
  }

  PolynomialFraction& operator-=(const value_type x)
  {
    dividend -= divisor*x;
    return *this;
  }

  PolynomialFraction& operator*=(const PolynomialFraction& p)
  {
    dividend *= p.dividend;
    divisor *= p.divisor;
    simplify();
    return *this;
  }

  PolynomialFraction& operator/=(const PolynomialFraction& p)
  {
    dividend *= p.divisor;
    divisor *= p.dividend;
    simplify();
    return *this;
  }

  PolynomialFraction& operator+=(const PolynomialFraction& p)
  {
    if (divisor != p.divisor)
    {
      dividend = dividend * p.divisor + p.dividend * divisor;
      divisor *= p.divisor;
    }
    else
    {
      dividend += p.dividend;
    }
    simplify();
    return *this;
  }

  PolynomialFraction& operator-=(const PolynomialFraction& p)
  {
    (*this) += -p;
    return *this;
  }

  PolynomialFraction operator-() const
  {
    PolynomialFraction result(*this);
    result *= -1.0;
    return result;
  }

  optimised_t optimise() const
  {
    return optimised_t(*this);
  }

  PolynomialFraction derivative(const variable_t& variable) const
  {
    const polynomial_t newDividend = dividend.derivative(variable)*divisor - dividend*divisor.derivative(variable);
    const polynomial_t newDivisor = divisor * divisor;
    return PolynomialFraction(newDividend, newDivisor);
  }
};

template<typename T>
std::ostream& operator<<(std::ostream& o, const PolynomialFraction<T>& f)
{
  o << "(" << f.getDividend() << ")/(" << f.getDivisor() << ")"; 
  return o;
}

}
#endif
