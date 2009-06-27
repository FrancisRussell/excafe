#include <utility>
#include <simple_cfd/numeric/polynomial.hpp>

namespace cfd
{

void Polynomial::addTerm(const double coefficient, const std::string& variable, const std::size_t exponent)
{
  addIndependentVariable(variable);
  coefficients.push_back(coefficient);

  for(exponent_map_t::iterator varIter(exponentMap.begin()); varIter!=exponentMap.end(); ++varIter)
  {
    varIter->second.push_back(varIter->first == variable ? exponent : 0);
  }
}

Polynomial::Polynomial()
{
}

Polynomial::Polynomial(const double constant)
{
  coefficients.push_back(constant);
}

Polynomial::Polynomial(const std::string& variable)
{
  addTerm(1.0, variable, 1);
}

Polynomial::Polynomial(const double coefficient, const std::string& variable)
{
  addTerm(coefficient, variable, 1);
}

Polynomial::Polynomial(const std::string& variable, const std::size_t exponent)
{
  addTerm(1.0, variable, exponent);
}

Polynomial::Polynomial(const double coefficient, const std::string& variable, const std::size_t exponent)
{
  addTerm(coefficient, variable, exponent);
}

void Polynomial::addIndependentVariable(const std::string& variable)
{
  if (exponentMap.find(variable) == exponentMap.end())
    exponentMap.insert(std::make_pair(variable, std::vector<std::size_t>(coefficients.size(), 0)));
}

}
