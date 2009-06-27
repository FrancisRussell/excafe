#ifndef SIMPLE_CFD_NUMERIC_POLYNOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_POLYNOMIAL_HPP

#include <vector>
#include <string>
#include <map>
#include <cstddef>

namespace cfd
{

class Polynomial
{
private:
  typedef std::map< std::string, std::vector<std::size_t> > exponent_map_t;
  std::vector<double> coefficients;
  exponent_map_t exponentMap;

  void multiply(const std::string& variable, const std::size_t exponent);
  void addTerm(const double coefficient, const std::string& variable, const std::size_t exponent);

public:
  Polynomial();

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

  Polynomial derivative(const std::size_t d) const;

  double operator()() const;
  double operator()(const double a) const;
  double operator()(const double a, const double b) const;
  double operator()(const double a, const double b, const double c) const;
};

}

#endif
