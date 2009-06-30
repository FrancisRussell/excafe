#include <utility>
#include <simple_cfd/numeric/polynomial.hpp>
#include <boost/lambda/lambda.hpp>
#include <cassert>
#include <ostream>

namespace cfd
{

Polynomial::Polynomial()
{
}

Polynomial::Polynomial(const Polynomial& p) : coefficients(p.coefficients), exponentMap(p.exponentMap)
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

void Polynomial::addTerm(const double coefficient, const std::string& variable, const std::size_t exponent)
{
  addIndependentVariable(variable);
  coefficients.push_back(coefficient);

  for(exponent_map_t::iterator varIter(exponentMap.begin()); varIter!=exponentMap.end(); ++varIter)
  {
    varIter->second.push_back(varIter->first == variable ? exponent : 0);
  }
}

void Polynomial::addConstant(const double constant)
{
  coefficients.push_back(constant);

  for(exponent_map_t::iterator varIter(exponentMap.begin()); varIter!=exponentMap.end(); ++varIter)
  {
    varIter->second.push_back(0);
  }
}

void Polynomial::addIndependentVariable(const std::string& variable)
{
  if (exponentMap.find(variable) == exponentMap.end())
    exponentMap.insert(std::make_pair(variable, std::vector<std::size_t>(coefficients.size(), 0)));
}

void Polynomial::addIndependentVariables(const Polynomial& p)
{
  for(exponent_map_t::const_iterator varIter(p.exponentMap.begin()); varIter!=p.exponentMap.end(); ++varIter)
  {
    addIndependentVariable(varIter->first);
  }
}

Polynomial& Polynomial::operator*=(const double x)
{
  transformCoefficients(boost::lambda::_1 * x);
  return *this;
}

Polynomial& Polynomial::operator/=(const double x)
{
  transformCoefficients(boost::lambda::_1 / x);
  return *this;
}

Polynomial& Polynomial::operator+=(const double x)
{
  addConstant(x);
  return *this;
}

Polynomial& Polynomial::operator-=(const double x)
{
  addConstant(-x);
  return *this;
}

Polynomial Polynomial::operator-() const
{
  Polynomial result(*this);
  result.transformCoefficients(-boost::lambda::_1);
  return result;
}

Polynomial& Polynomial::operator*=(const Polynomial& p)
{
  Polynomial result;

  for(std::size_t pTerm=0; pTerm<p.numTerms(); ++pTerm)
  {
    Polynomial multiplied(*this);
    multiplied.transformCoefficients(boost::lambda::_1 * p.coefficients[pTerm]);
    
    for(exponent_map_t::const_iterator pExpIter(p.exponentMap.begin()); pExpIter!=p.exponentMap.end(); ++pExpIter)
    {
      multiplied.addIndependentVariable(pExpIter->first);
      multiplied.transformExponents(pExpIter->first, boost::lambda::_1 + pExpIter->second[pTerm]);
    }

    result += multiplied;
  }

  result.swap(*this);

  return *this;
}

Polynomial& Polynomial::operator+=(const Polynomial& p)
{
  addIndependentVariables(p);

  coefficients.insert(coefficients.end(), p.coefficients.begin(), p.coefficients.end());

  for(exponent_map_t::iterator varIter(exponentMap.begin()); varIter!=exponentMap.end(); ++varIter)
  {
    const exponent_map_t::const_iterator pIter = p.exponentMap.find(varIter->first);
    if (pIter != p.exponentMap.end())
    {
      varIter->second.insert(varIter->second.end(), pIter->second.begin(), pIter->second.end());
    }
    else
    {
      varIter->second.resize(numTerms(), 0.0);
    }
  }

  return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& p)
{
  *this += -p;
  return *this;
}

std::size_t Polynomial::numTerms() const
{
  return coefficients.size();
}

void Polynomial::swap(Polynomial& p)
{
  coefficients.swap(p.coefficients);
  exponentMap.swap(p.exponentMap);
}

std::ostream& operator<<(std::ostream& out, const Polynomial& p)
{
  for(std::size_t pTerm=0; pTerm<p.numTerms(); ++pTerm)
  {
    out << p.coefficients[pTerm];

    for(Polynomial::exponent_map_t::const_iterator eIter(p.exponentMap.begin()); eIter!=p.exponentMap.end(); ++eIter)
    {
      if (eIter->second[pTerm] > 0)
      {
        out << " * ";
        out << eIter->first;
        
        if (eIter->second[pTerm] != 1)
          out << "^" << eIter->second[pTerm];
      }
    }

    if (pTerm < p.numTerms() - 1)
      out << " + ";
  }

  return out;
}

}
