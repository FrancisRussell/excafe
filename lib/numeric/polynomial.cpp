#include <utility>
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/numeric/monomial.hpp>
#include <simple_cfd/numeric/optimised_polynomial.hpp>
#include <boost/lambda/lambda.hpp>
#include <cassert>
#include <algorithm>
#include <ostream>

namespace cfd
{

Polynomial::Polynomial()
{
}

Polynomial::Polynomial(const Polynomial& p) : independentVariables(p.independentVariables), coefficients(p.coefficients)
{
}

Polynomial::Polynomial(const double constant)
{
  addConstant(constant);
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
  coefficients[Monomial(variable, exponent)] += coefficient;
  cleanZeros();
}

void Polynomial::addConstant(const double constant)
{
  coefficients[Monomial()] += constant;
  cleanZeros();
}

void Polynomial::addIndependentVariable(const std::string& variable)
{
  independentVariables.insert(variable);
}

void Polynomial::addIndependentVariables(const Polynomial& p)
{
  independentVariables.insert(p.independentVariables.begin(), p.independentVariables.end());
}

void Polynomial::addIndependentVariables(const Monomial& m)
{
  const std::set<std::string> mVariables(m.getVariables()); 
  independentVariables.insert(mVariables.begin(), mVariables.end());
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

Polynomial& Polynomial::operator*=(const Monomial& m)
{
  addIndependentVariables(m);
  std::map<Monomial, double> newCoefficients;

  for(coefficient_map_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
    newCoefficients.insert(std::make_pair(cIter->first * m, cIter->second));

  coefficients.swap(newCoefficients);
  return *this;
}

Polynomial Polynomial::operator*(const Monomial& m) const
{
  Polynomial result(*this);
  result *= m;
  return result;
}

Polynomial& Polynomial::operator*=(const Polynomial& p)
{
  Polynomial result;

  for(coefficient_map_t::const_iterator cIter(p.coefficients.begin()); cIter!=p.coefficients.end(); ++cIter)
  {
    result += (*this) * cIter->first * cIter->second;
  }

  result.swap(*this);
  return *this;
}

Polynomial& Polynomial::operator+=(const Polynomial& p)
{
  addIndependentVariables(p);

  for(coefficient_map_t::const_iterator cIter(p.coefficients.begin()); cIter!=p.coefficients.end(); ++cIter)
    coefficients[cIter->first] += cIter->second;

  cleanZeros();
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

std::set<std::string> Polynomial::getIndependentVariables() const
{
  return independentVariables;
}

void Polynomial::swap(Polynomial& p)
{
  coefficients.swap(p.coefficients);
  independentVariables.swap(p.independentVariables);
}

void Polynomial::cleanZeros()
{
  std::set<Monomial> zeroMonomials;

  for(Polynomial::coefficient_map_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
  {
    if (cIter->second == 0.0)
      zeroMonomials.insert(cIter->first);
  }

  for(std::set<Monomial>::const_iterator mIter(zeroMonomials.begin()); mIter!=zeroMonomials.end(); ++mIter)
    coefficients.erase(*mIter);
}

Polynomial::iterator Polynomial::begin()
{
  return coefficients.begin();
}

Polynomial::iterator Polynomial::end()
{
  return coefficients.end();
}

Polynomial::const_iterator Polynomial::begin() const
{
  return coefficients.begin();
}

Polynomial::const_iterator Polynomial::end() const
{
  return coefficients.end();
}

Polynomial Polynomial::derivative(const std::string& variable) const
{
  Polynomial result;
  result.addIndependentVariables(*this);

  for(Polynomial::coefficient_map_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
  {
    const std::pair<double, Monomial> mDerivative(cIter->first.derivative(variable));
    result.coefficients.insert(std::make_pair(mDerivative.second, mDerivative.first * cIter->second));
  }

  result.cleanZeros();
  return result;
}

OptimisedPolynomial Polynomial::optimise() const
{
  return OptimisedPolynomial(*this);
}

void Polynomial::checkConsistent() const
{
  // Checks that independentVariables is a superset of all variables referenced
  // in monomials.
  std::set<std::string> referencedVariables;

  for(Polynomial::coefficient_map_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
  {
    const std::set<std::string> monomialVars(cIter->first.getVariables());
    referencedVariables.insert(monomialVars.begin(), monomialVars.end()); 
  }

  assert(std::includes(independentVariables.begin(), independentVariables.end(),
         referencedVariables.begin(), referencedVariables.end()));
}

std::ostream& operator<<(std::ostream& out, const Polynomial& p)
{
  for(Polynomial::coefficient_map_t::const_iterator cIter(p.coefficients.begin()); cIter!=p.coefficients.end(); ++cIter)
  {
    if (cIter != p.coefficients.begin())
      out << " + ";

    const bool renderCoefficient = cIter->second != 1.0 || cIter->first.isOne();
    const bool renderMonomial = !cIter->first.isOne();

    if (renderCoefficient)
      out << cIter->second; 
      
    if (renderCoefficient && renderMonomial)
      out << "*";
      
    if (renderMonomial) 
      out << cIter->first;
  }
  
  if (p.coefficients.empty())
    out << 0.0;

  return out;
}

}
