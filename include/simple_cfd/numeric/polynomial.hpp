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
#include <boost/utility.hpp>
#include <simple_cfd_fwd.hpp>
#include <numeric/monomial.hpp>
#include <numeric/optimised_polynomial.hpp>
#include <util/lazy_copy.hpp>

namespace cfd
{

template<typename V>
class Polynomial : boost::addable<Polynomial<V>, double,
                   boost::subtractable<Polynomial<V>, double,
                   boost::dividable<Polynomial<V>, double,
                   boost::multipliable<Polynomial<V>, double, 
                   boost::addable<Polynomial<V>,
                   boost::subtractable< Polynomial<V>,
                   boost::multipliable< Polynomial<V>,
                   boost::equality_comparable< Polynomial<V>
                   > > > > > > > >
{
public:
  typedef V variable_t;
  typedef OptimisedPolynomial<variable_t> optimised_t;

private:
  typedef std::map<Monomial<variable_t>, double> coefficient_map_t;
  util::LazyCopy<coefficient_map_t> coefficients;

  void addTerm(const double coefficient, const variable_t& variable, const std::size_t exponent)
  {
    (*coefficients)[Monomial<variable_t>(variable, exponent)] += coefficient;
  }

  void addConstant(const double constant)
  {
    (*coefficients)[Monomial<variable_t>()] += constant;
  }

  // We currently only call this from public methods
  void cleanZeros()
  {
    typedef typename coefficient_map_t::iterator coeff_map_iter;

    // We need this construction because we invalidate the current iterator on erase
    coeff_map_iter currentIter = coefficients->begin();
    while(currentIter != coefficients->end())
    {
      const coeff_map_iter nextIter = boost::next(currentIter);

      if (currentIter->second == 0.0)
        coefficients->erase(currentIter);

      currentIter = nextIter;
    }
  }

  void addMonomial(const double coefficient, const Monomial<variable_t>& m)
  {
    (*coefficients)[m] += coefficient;
  }

  template<typename UnaryFunction>
  void transformCoefficients(const UnaryFunction& f)
  {
    for(typename coefficient_map_t::iterator cIter(coefficients->begin()); cIter!=coefficients->end(); ++cIter)
      cIter->second = f(cIter->second);
  }

public:
  typedef typename coefficient_map_t::iterator iterator;
  typedef typename coefficient_map_t::const_iterator const_iterator;

  Polynomial()
  {
  }

  Polynomial(const Polynomial& p) : coefficients(p.coefficients)
  {
  }

  Polynomial(const double constant)
  {
    addConstant(constant);
    cleanZeros();
  }

  Polynomial(const variable_t& variable)
  {
    addTerm(1.0, variable, 1);
    cleanZeros();
  }

  Polynomial(const double coefficient, const variable_t& variable)
  {
    addTerm(coefficient, variable, 1);
    cleanZeros();
  }

  Polynomial(const variable_t& variable, const std::size_t exponent)
  {
    addTerm(1.0, variable, exponent);
    cleanZeros();
  }

  Polynomial(const double coefficient, const variable_t& variable, const std::size_t exponent)
  {
    addTerm(coefficient, variable, exponent);
    cleanZeros();
  }

  iterator begin()
  {
    return coefficients->begin();
  }

  iterator end()
  {
    return coefficients->end();
  }

  const_iterator begin() const
  {
    return coefficients->begin();
  }

  const_iterator end() const
  {
    return coefficients->end();
  }

  void replaceVariable(const variable_t& from, const variable_t& to)
  {
    coefficient_map_t oldCoefficients;
    std::swap(coefficients, oldCoefficients);

    BOOST_FOREACH(const typename coefficient_map_t::value_type& monomialMapping, oldCoefficients)
    {
      // We specifically use addTerm here, to avoid issues when previously distinct monomials
      // become identical after variable substitution
      addTerm(monomialMapping.first.substitute(from, to), monomialMapping.second);
    }
    cleanZeros();
  }

  std::set<variable_t> getVariables() const
  {
    std::set<variable_t> result;

    for(typename coefficient_map_t::const_iterator cIter(coefficients->begin()); cIter!=coefficients->end(); ++cIter)
    {
      const std::set<variable_t> monomialVariables = cIter->first.getVariables();
      result.insert(monomialVariables.begin(), monomialVariables.end());
    }

    return result;
  }

  void checkConsistent() const
  {
    // NOOP
  }

  bool operator==(const Polynomial& p) const
  {
    return coefficients == p.coefficients;
  }

  Polynomial& operator*=(const double x)
  {
    transformCoefficients(boost::lambda::_1 * x);
    cleanZeros();
    return *this;
  }

  Polynomial& operator/=(const double x)
  {
    transformCoefficients(boost::lambda::_1 / x);
    cleanZeros();
    return *this;
  }

  Polynomial& operator+=(const double x)
  {
    addConstant(x);
    cleanZeros();
    return *this;
  }

  Polynomial& operator-=(const double x)
  {
    addConstant(-x);
    cleanZeros();
    return *this;
  }

  Polynomial& operator*=(const Polynomial& b)
  {
    Polynomial a;
    std::swap(a, *this);
  
    BOOST_FOREACH(const typename coefficient_map_t::value_type& aMapping, a.coefficients.cref())
    {
      BOOST_FOREACH(const typename coefficient_map_t::value_type& bMapping, b.coefficients.cref())
      {
        addMonomial(aMapping.second * bMapping.second, aMapping.first * bMapping.first);
      }
    }

    cleanZeros();
    return *this;
  }

  Polynomial& operator+=(const Polynomial& p)
  {
    typename coefficient_map_t::iterator coeffIter = coefficients->begin();

    BOOST_FOREACH(const typename coefficient_map_t::value_type& pMapping, *p.coefficients)
    {
      while(coeffIter != coefficients->end() && coeffIter->first < pMapping.first)
        ++coeffIter;

      if (coeffIter == coefficients->end() || !(coeffIter->first == pMapping.first))
        coefficients->insert(coeffIter, pMapping);
      else
        coeffIter->second += pMapping.second;
    }

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
    result *= -1.0;
    return result;
  }

  Polynomial derivative(const variable_t& variable) const
  {
    Polynomial result;
    for(typename coefficient_map_t::const_iterator cIter(coefficients->begin()); cIter!=coefficients->end(); ++cIter)
    {
      const std::pair< double, Monomial<variable_t> > mDerivative(cIter->first.derivative(variable));
      result.addMonomial(cIter->second * mDerivative.first, mDerivative.second);
    }
  
    result.cleanZeros();
    return result;
  }

  Polynomial substituteValue(const variable_t& variable, const double value) const
  {
    Polynomial result;
  
    for(typename coefficient_map_t::const_iterator cIter(coefficients->begin()); cIter!=coefficients->end(); ++cIter)
    {
      const std::pair< double, Monomial<variable_t> > mBound(cIter->first.substituteValue(variable, value));
      result.addMonomial(cIter->second * mBound.first, mBound.second);
    }
  
    result.cleanZeros();
    return result;
  }

  optimised_t optimise() const
  {
    return optimised_t(*this);
  }

  std::size_t numTerms() const
  {
    return coefficients->size();
  }

  std::ostream& write(std::ostream& out) const
  {
    for(typename coefficient_map_t::const_iterator cIter(coefficients->begin()); cIter!=coefficients->end(); ++cIter)
    {
      if (cIter != coefficients->begin())
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
    
    if (coefficients->empty())
      out << 0.0;
  
    return out;
  }

  void swap(Polynomial& p)
  {
    coefficients.swap(p.coefficients);
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
