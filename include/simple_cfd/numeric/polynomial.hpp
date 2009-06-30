#ifndef SIMPLE_CFD_NUMERIC_POLYNOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_POLYNOMIAL_HPP

#include <vector>
#include <string>
#include <map>
#include <cstddef>
#include <algorithm>
#include <cassert>
#include <iosfwd>
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

  typedef std::map< std::string, std::vector<std::size_t> > exponent_map_t;
  std::vector<double> coefficients;
  exponent_map_t exponentMap;

  void multiply(const std::string& variable, const std::size_t exponent);
  void addTerm(const double coefficient, const std::string& variable, const std::size_t exponent);
  void addConstant(const double constant);
  void addIndependentVariables(const Polynomial& p);

  template<typename UnaryFunction>
  void transformCoefficients(const UnaryFunction& f)
  {
    std::transform(coefficients.begin(), coefficients.end(), coefficients.begin(), f);
  }

  template<typename UnaryFunction>
  void transformExponents(const std::string& variable, const UnaryFunction& f)
  {
    exponent_map_t::iterator eIter(exponentMap.find(variable));
    assert(eIter != exponentMap.end());
    std::transform(eIter->second.begin(), eIter->second.end(), eIter->second.begin(), f);
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

  Polynomial& operator*=(const double x);
  Polynomial& operator/=(const double x);
  Polynomial& operator+=(const double x);
  Polynomial& operator-=(const double x);

  Polynomial& operator*=(const Polynomial& p);
  Polynomial& operator+=(const Polynomial& p);
  Polynomial& operator-=(const Polynomial& p);

  Polynomial operator-() const;

  Polynomial derivative(const std::size_t d) const;
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
