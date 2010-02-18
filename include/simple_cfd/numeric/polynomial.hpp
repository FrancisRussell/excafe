#ifndef SIMPLE_CFD_NUMERIC_POLYNOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_POLYNOMIAL_HPP

#include <vector>
#include <map>
#include <set>
#include <cstddef>
#include <algorithm>
#include <cassert>
#include <utility>
#include <iosfwd>
#include <boost/operators.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/foreach.hpp>
#include <simple_cfd_fwd.hpp>
#include <numeric/monomial.hpp>
#include <numeric/optimised_polynomial.hpp>

namespace cfd
{

template<typename V>
class Polynomial : boost::addable2<Polynomial<V>, double,
                   boost::subtractable2<Polynomial<V>, double,
                   boost::dividable2<Polynomial<V>, double,
                   boost::multipliable2<Polynomial<V>, double, 
                   boost::addable< Polynomial<V>,
                   boost::subtractable< Polynomial<V>,
                   boost::multipliable< Polynomial<V>
                   > > > > > > >
{
public:
  typedef V variable_t;

private:
  typedef std::map<Monomial<variable_t>, double> coefficient_map_t;
  std::set<variable_t> independentVariables;
  coefficient_map_t coefficients;

  void addTerm(const double coefficient, const variable_t& variable, const std::size_t exponent)
  {
    addIndependentVariable(variable);
    coefficients[Monomial<variable_t>(variable, exponent)] += coefficient;
    cleanZeros();
  }

  void addTerm(const Monomial<variable_t>& m, const double coefficient)
  {
    addIndependentVariables(m);
    coefficients[m] += coefficient;
  }

  void addConstant(const double constant)
  {
    coefficients[Monomial<variable_t>()] += constant;
    cleanZeros();
  }

  void addIndependentVariables(const Polynomial& p)
  {
    independentVariables.insert(p.independentVariables.begin(), p.independentVariables.end());
  }

  void addIndependentVariables(const Monomial<variable_t>& m)
  {
    const std::set<variable_t> mVariables(m.getVariables()); 
    independentVariables.insert(mVariables.begin(), mVariables.end());
  }

  void cleanZeros()
  {
    std::set< Monomial<variable_t> > zeroMonomials;
  
    for(typename coefficient_map_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
    {
      if (cIter->second == 0.0)
        zeroMonomials.insert(cIter->first);
    }
  
    for(typename std::set< Monomial<variable_t> >::const_iterator mIter(zeroMonomials.begin()); mIter!=zeroMonomials.end(); ++mIter)
      coefficients.erase(*mIter);
  }


  void addMonomial(const double coefficient, const Monomial<variable_t>& m)
  {
    coefficients[m] += coefficient;
  }

  Polynomial operator*(const Monomial<variable_t>& m) const
  {
    Polynomial result(*this);
    result *= m;
    return result;
  }

  Polynomial& operator*=(const Monomial<variable_t>& m)
  {
    addIndependentVariables(m);
    std::map<Monomial<variable_t>, double> newCoefficients;
  
    for(typename coefficient_map_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
      newCoefficients.insert(std::make_pair(cIter->first * m, cIter->second));
  
    coefficients.swap(newCoefficients);
    return *this;
  }


  template<typename UnaryFunction>
  void transformCoefficients(const UnaryFunction& f)
  {
    for(typename coefficient_map_t::iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
      cIter->second = f(cIter->second);

    cleanZeros();
  }

public:
  typedef typename coefficient_map_t::iterator iterator;
  typedef typename coefficient_map_t::const_iterator const_iterator;

  Polynomial()
  {
  }

  Polynomial(const Polynomial& p) : independentVariables(p.independentVariables), coefficients(p.coefficients)
  {
  }

  explicit Polynomial(const double constant)
  {
    addConstant(constant);
  }

  explicit Polynomial(const variable_t& variable)
  {
    addTerm(1.0, variable, 1);
  }

  explicit Polynomial(const double coefficient, const variable_t& variable)
  {
    addTerm(coefficient, variable, 1);
  }

  explicit Polynomial(const variable_t& variable, const std::size_t exponent)
  {
    addTerm(1.0, variable, exponent);
  }

  explicit Polynomial(const double coefficient, const variable_t& variable, const std::size_t exponent)
  {
    addTerm(coefficient, variable, exponent);
  }

  iterator begin()
  {
    return coefficients.begin();
  }

  iterator end()
  {
    return coefficients.end();
  }

  const_iterator begin() const
  {
    return coefficients.begin();
  }

  const_iterator end() const
  {
    return coefficients.end();
  }

  void addIndependentVariable(const variable_t& variable)
  {
    independentVariables.insert(variable);
  }

  void replaceIndependentVariable(const variable_t& from, const variable_t& to)
  {
    const typename std::set<variable_t>::const_iterator fromIter = independentVariables.find(from);

    if (fromIter != independentVariables.end())
    {
      independentVariables.erase(fromIter);
      independentVariables.insert(to);

      coefficient_map_t oldCoefficients;
      std::swap(coefficients, oldCoefficients);

      BOOST_FOREACH(typename coefficient_map_t::value_type monomialMapping, oldCoefficients)
      {
        // We specifically use addTerm here, to avoid issues when previously distinct monomials
        // become identical after variable substitution
        addTerm(monomialMapping.first.substitute(from, to), monomialMapping.second);
      }
    }
  }

  std::set<variable_t> getIndependentVariables() const
  {
    return independentVariables;
  }

  void checkConsistent() const
  {
    // Checks that independentVariables is a superset of all variables referenced
    // in monomials.
    std::set<variable_t> referencedVariables;
  
    for(typename coefficient_map_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
    {
      const std::set<variable_t> monomialVars(cIter->first.getVariables());
      referencedVariables.insert(monomialVars.begin(), monomialVars.end()); 
    }
  
    assert(std::includes(independentVariables.begin(), independentVariables.end(),
           referencedVariables.begin(), referencedVariables.end()));
  }

  Polynomial& operator*=(const double x)
  {
    transformCoefficients(boost::lambda::_1 * x);
    return *this;
  }

  Polynomial& operator/=(const double x)
  {
    transformCoefficients(boost::lambda::_1 / x);
    return *this;
  }

  Polynomial& operator+=(const double x)
  {
    addConstant(x);
    return *this;
  }

  Polynomial& operator-=(const double x)
  {
    addConstant(-x);
    return *this;
  }

  Polynomial& operator*=(const Polynomial& p)
  {
    Polynomial result;
  
    for(typename coefficient_map_t::const_iterator cIter(p.coefficients.begin()); cIter!=p.coefficients.end(); ++cIter)
    {
      result += (*this) * cIter->first * cIter->second;
    }
  
    result.swap(*this);
    return *this;
  }

  Polynomial& operator+=(const Polynomial& p)
  {
    addIndependentVariables(p);
  
    for(typename coefficient_map_t::const_iterator cIter(p.coefficients.begin()); cIter!=p.coefficients.end(); ++cIter)
      coefficients[cIter->first] += cIter->second;
  
    cleanZeros();
    return *this;
  }

  Polynomial& operator-=(const Polynomial& p)
  {
    *this += -p;
    return *this;
  }

  Polynomial operator-() const
  {
    Polynomial result(*this);
    result.transformCoefficients(-boost::lambda::_1);
    return result;
  }


  Polynomial derivative(const variable_t& variable) const
  {
    Polynomial result;
    result.addIndependentVariables(*this);
  
    for(typename coefficient_map_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
    {
      const std::pair< double, Monomial<variable_t> > mDerivative(cIter->first.derivative(variable));
      result.addMonomial(cIter->second * mDerivative.first, mDerivative.second);
    }
  
    result.cleanZeros();
    return result;
  }

  OptimisedPolynomial<variable_t> optimise() const
  {
    return OptimisedPolynomial<variable_t>(*this);
  }

  std::size_t numTerms() const
  {
    return coefficients.size();
  }

  std::ostream& write(std::ostream& out) const
  {
    for(typename coefficient_map_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
    {
      if (cIter != coefficients.begin())
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
    
    if (coefficients.empty())
      out << 0.0;
  
    return out;
  }

  void swap(Polynomial& p)
  {
    coefficients.swap(p.coefficients);
    independentVariables.swap(p.independentVariables);
  }
};

template<typename V>
std::ostream& operator<<(std::ostream& out, const Polynomial<V>& p)
{
  return p.write(out);
}

}

namespace std
{
  template<typename V>
  void swap(cfd::Polynomial<V>& a, cfd::Polynomial<V>& b)
  {
    a.swap(b);
  }
}

#endif
