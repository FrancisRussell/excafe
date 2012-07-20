#ifndef EXCAFE_NUMERIC_POLYNOMIAL_FRACTION_HPP
#define EXCAFE_NUMERIC_POLYNOMIAL_FRACTION_HPP

#include <cstddef>
#include <cstdlib>
#include <map>
#include <boost/operators.hpp>
#include "polynomial.hpp"
#include "numeric_fwd.hpp"
#include "expression.hpp"
#include "expression_visitor.hpp"
#include "traits.hpp"
#include <excafe/exception.hpp>
#include <excafe/mp/float.hpp>
#include <excafe/mp/rational.hpp>
#include <excafe/mp/integer.hpp>
#include <ostream>

namespace excafe
{

template<typename V>
class PolynomialFraction : public NumericExpression<V>,
                           boost::arithmetic<PolynomialFraction<V>, typename Polynomial<V>::value_type,
                           boost::arithmetic<PolynomialFraction<V>,
                           boost::totally_ordered< PolynomialFraction<V>
                           > > >
{
public:
  static const bool supports_abs = false;

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

  PolynomialFraction(const mp::Float& constant) : dividend(constant), divisor(1.0)
  {
  }

  PolynomialFraction(const mp::Integer& constant) : dividend(constant), divisor(1.0)
  {
  }

  PolynomialFraction(const mp::Rational& r) : dividend(r.getNumerator()), divisor(r.getDenominator())
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

  void accept(NumericExpressionVisitor<variable_t>& v) const
  {
    dividend.accept(v);

    if (divisor != polynomial_t(1.0))
    {
      divisor.accept(v);
      v.visitExponent(-1);
      v.postProduct(2);
    }
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

  PolynomialFraction reciprocal() const
  {
    assert(dividend != 0);
    return PolynomialFraction(divisor, dividend);
  }

  PolynomialFraction pow(const int n) const
  {
    const int absExponent = std::abs(n);
    const PolynomialFraction raised(dividend.pow(absExponent), divisor.pow(absExponent));

    if (n >= 0)
      return raised;
    else
      return raised.reciprocal();
  }

  optimised_t optimise() const
  {
    return optimised_t(*this);
  }

  std::size_t degree(const variable_t& variable) const
  {
    if (divisor.degree(variable)>0)
      CFD_EXCEPTION("Cannot compute degree of variable in numerator and denominator of fraction.");
    else
      return dividend.degree(variable);
  }

  PolynomialFraction derivative(const variable_t& variable) const
  {
    const polynomial_t newDividend = dividend.derivative(variable)*divisor - dividend*divisor.derivative(variable);
    const polynomial_t newDivisor = divisor * divisor;
    return PolynomialFraction(newDividend, newDivisor);
  }

  void write(std::ostream& o) const
  {
    const bool noDivisor = (divisor == polynomial_t(1.0));

    if (!noDivisor)
      o << "(";

    o << dividend;

    if (!noDivisor)
      o << ") * (" << divisor << ")^-1";
  }
};

template<typename V>
PolynomialFraction<V> pow(const PolynomialFraction<V>& f, const int n)
{
  return f.pow(n);
}

template<typename T>
std::ostream& operator<<(std::ostream& o, const PolynomialFraction<T>& f)
{
  f.write(o);
  return o;
}

}
#endif
