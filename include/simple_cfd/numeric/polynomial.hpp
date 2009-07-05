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
#include <boost/operators.hpp>
#include <simple_cfd_fwd.hpp>
#include <numeric/monomial.hpp>

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

  typedef std::map<Monomial, double> coefficient_map_t;
  std::set<std::string> independentVariables;
  coefficient_map_t coefficients;

  void multiply(const std::string& variable, const std::size_t exponent);
  void addTerm(const double coefficient, const std::string& variable, const std::size_t exponent);
  void addConstant(const double constant);
  void addIndependentVariables(const Polynomial& p);
  void addIndependentVariables(const Monomial& m);
  void cleanZeros();

  Polynomial& operator*=(const Monomial& m);
  Polynomial operator*(const Monomial& m) const;

  template<typename UnaryFunction>
  void transformCoefficients(const UnaryFunction& f)
  {
    for(coefficient_map_t::iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
      cIter->second = f(cIter->second);

    cleanZeros();
  }

public:
  typedef coefficient_map_t::iterator iterator;
  typedef coefficient_map_t::const_iterator const_iterator;

  Polynomial();
  Polynomial(const Polynomial& p);

  explicit Polynomial(const double constant);
  explicit Polynomial(const std::string& variable);
  explicit Polynomial(const double coefficient, const std::string& variable);
  explicit Polynomial(const std::string& variable, const std::size_t exponent);
  explicit Polynomial(const double coefficient, const std::string& variable, const std::size_t exponent);

  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;

  void addIndependentVariable(const std::string& s);
  std::set<std::string> getIndependentVariables() const;
  void checkConsistent() const;

  Polynomial& operator*=(const double x);
  Polynomial& operator/=(const double x);
  Polynomial& operator+=(const double x);
  Polynomial& operator-=(const double x);

  Polynomial& operator*=(const Polynomial& p);
  Polynomial& operator+=(const Polynomial& p);
  Polynomial& operator-=(const Polynomial& p);

  Polynomial operator-() const;

  Polynomial derivative(const std::string& variable) const;
  OptimisedPolynomial optimise() const;
  std::size_t numTerms() const;
  void swap(Polynomial& p);
};

std::ostream& operator<<(std::ostream& o, const Polynomial& p);

}

#endif
