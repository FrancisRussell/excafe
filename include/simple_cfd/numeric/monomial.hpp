#ifndef SIMPLE_CFD_NUMERIC_MONOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_MONOMIAL_HPP

#include <map>
#include <string>
#include <cstddef>
#include <set>
#include <iosfwd>
#include <utility>

namespace cfd
{

class Monomial
{
private:
  friend std::ostream& operator<<(std::ostream& o, const Monomial& m);
  std::map<std::string, std::size_t> exponents;

public:
  Monomial();
  Monomial(const Monomial& m);
  Monomial(const std::string& variable, const std::size_t exponent);
  Monomial& operator*=(const Monomial& m);
  Monomial operator*(const Monomial& m) const;
  bool operator==(const Monomial& m) const;
  bool operator<(const Monomial& m) const;
  Monomial& operator=(const Monomial& m);

  bool isOne() const;
  std::set<std::string> getVariables() const;
  std::size_t getExponent(const std::string& variable) const;
  std::pair<double, Monomial> derivative(const std::string& variable) const;
};

std::ostream& operator<<(std::ostream& o, const Monomial& p);

}

#endif
