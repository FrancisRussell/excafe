#ifndef SIMPLE_CFD_NUMERIC_POLYNOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_POLYNOMIAL_HPP

#include <vector>
#include <string>
#include <map>
#include <set>
#include <cstddef>
#include <algorithm>
#include <cassert>
#include <iosfwd>
#include <numeric/monomial.hpp>
#include <boost/operators.hpp>

namespace cfd
{

class Polynomial : boost::addable2<Polynomial, double,
                   boost::subtractable2<Polynomial, double,
                   boost::dividable2<Polynomial, double,
                   boost::multipliable2<Polynomial, double, 
                   boost::addable< Polynomial,
                   boost::subtractable< Polynomial,
                   boost::multipliable< Polynomial
                   > > > > > > >
{
private:
  friend std::ostream& operator<<(std::ostream& o, const Polynomial& p);

  typedef std::map<detail::Monomial, double> coefficient_map_t;
  std::set<std::string> independentVariables;
  coefficient_map_t coefficients;

  void multiply(const std::string& variable, const std::size_t exponent);
  void addTerm(const double coefficient, const std::string& variable, const std::size_t exponent);
  void addConstant(const double constant);
  void addIndependentVariables(const Polynomial& p);
  void cleanZeros();

  Polynomial& operator*=(const detail::Monomial& m);
  Polynomial operator*(const detail::Monomial& m) const;

  template<typename UnaryFunction>
  void transformCoefficients(const UnaryFunction& f)
  {
    for(coefficient_map_t::iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
      cIter->second = f(cIter->second);

    cleanZeros();
  }

public:
  Polynomial();
  Polynomial(const Polynomial& p);

  explicit Polynomial(const double constant);
  explicit Polynomial(const std::string& variable);
  explicit Polynomial(const double coefficient, const std::string& variable);
  explicit Polynomial(const std::string& variable, const std::size_t exponent);
  explicit Polynomial(const double coefficient, const std::string& variable, const std::size_t exponent);

  void addIndependentVariable(const std::string& s);
  std::set<std::string> getIndependentVariables() const;

  Polynomial& operator*=(const double x);
  Polynomial& operator/=(const double x);
  Polynomial& operator+=(const double x);
  Polynomial& operator-=(const double x);

  Polynomial& operator*=(const Polynomial& p);
  Polynomial& operator+=(const Polynomial& p);
  Polynomial& operator-=(const Polynomial& p);

  Polynomial operator-() const;

  Polynomial derivative(const std::string& variable) const;
  std::size_t numTerms() const;
  void swap(Polynomial& p);

  double operator()() const;
  double operator()(const double a) const;
  double operator()(const double a, const double b) const;
  double operator()(const double a, const double b, const double c) const;
};

std::ostream& operator<<(std::ostream& o, const Polynomial& p);

}

#endif
