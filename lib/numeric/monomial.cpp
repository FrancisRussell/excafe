#include <numeric/monomial.hpp>
#include <ostream>
#include <cstddef>
#include <map>
#include <string>
#include <utility>

namespace cfd
{

Monomial::Monomial()
{
}

Monomial::Monomial(const Monomial& m) : exponents(m.exponents)
{
}

Monomial::Monomial(const std::string& variable, const std::size_t exponent)
{
  exponents.insert(std::make_pair(variable, exponent));
}

Monomial& Monomial::operator*=(const Monomial& m)
{
  for (std::map<std::string, std::size_t>::const_iterator expIter(m.exponents.begin()); expIter!=m.exponents.end(); ++expIter)
    exponents[expIter->first] += expIter->second;

  return *this;
}

Monomial Monomial::operator*(const Monomial& m) const
{
  Monomial result(*this);
  result *= m;
  return result;
}

bool Monomial::operator==(const Monomial& m) const
{
  return exponents == m.exponents;
}

bool Monomial::operator<(const Monomial& m) const
{
  return exponents < m.exponents;
}

Monomial& Monomial::operator=(const Monomial& m)
{
  exponents = m.exponents;
  return *this;
}

std::set<std::string> Monomial::getVariables() const
{
  std::set<std::string> variables;

  for (std::map<std::string, std::size_t>::const_iterator expIter(exponents.begin()); expIter!=exponents.end(); ++expIter)
    variables.insert(expIter->first);

  return variables;
}

bool Monomial::isOne() const
{
  return exponents.empty();
}

std::size_t Monomial::getExponent(const std::string& variable) const
{
  const std::map<std::string, std::size_t>::const_iterator varIter(exponents.find(variable));

  if (varIter == exponents.end())
  {
    return 0;
  }
  else
  {
    return varIter->second;
  }
}

std::pair<double, Monomial> Monomial::derivative(const std::string& variable) const
{
  const double coefficient = getExponent(variable);
  Monomial result;

  for (std::map<std::string, std::size_t>::const_iterator eIter(exponents.begin()); eIter!=exponents.end(); ++eIter)
  {
    if (eIter->first != variable)
    {
      result.exponents.insert(*eIter);
    }
    else if (eIter->second > 1)
    {
      result.exponents.insert(std::make_pair(eIter->first, eIter->second - 1));
    }
  }

  return std::make_pair(coefficient, result);
}

void Monomial::swap(Monomial& m)
{
  exponents.swap(m.exponents); 
}

std::ostream& operator<<(std::ostream& out, const Monomial& m)
{
  for (std::map<std::string, std::size_t>::const_iterator eIter(m.exponents.begin()); eIter!=m.exponents.end(); ++eIter)
  {
    out << eIter->first;

    if (eIter->second != 1)
      out << "^{" << eIter->second << "}";
  }
  
  if (m.exponents.empty())
    out << "1.0";

  return out;
}

}

namespace std
{
  template<>
  void swap(cfd::Monomial& a, cfd::Monomial& b)
  {
    a.swap(b);
  }
}
