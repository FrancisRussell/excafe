#include <numeric/monomial.hpp>
#include <ostream>

namespace cfd
{

namespace detail
{

Monomial::Monomial()
{
}

Monomial::Monomial(const Monomial& m) : exponents(m.exponents)
{
}

Monomial::Monomial(const std::string& variable, const std::size_t exponent)
{
  exponents[variable] = exponent;
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

}
