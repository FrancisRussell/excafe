#ifndef SIMPLE_CFD_NUMERIC_POLYNOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_POLYNOMIAL_HPP

#include <vector>
#include <map>
#include <set>
#include <cstddef>
#include <algorithm>
#include <cassert>
#include <utility>
#include <numeric>
#include <iosfwd>
#include <boost/operators.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <simple_cfd_fwd.hpp>
#include "monomial.hpp"
#include "optimised_polynomial.hpp"
#include "cln_wrapper.hpp"
#include "value_map.hpp"
#include "expression.hpp"
#include "expression_visitor.hpp"
#include <simple_cfd/util/lazy_copy.hpp>
#include <simple_cfd/exception.hpp>

namespace cfd
{

template<typename V>
class Polynomial : public NumericExpression<V>,
                   boost::arithmetic<Polynomial<V>, double,
                   boost::additive<Polynomial<V>,
                   boost::multipliable< Polynomial<V>,
                   boost::totally_ordered< Polynomial<V>
                   > > > >
{
private:
  static const std::size_t precision = 32;

public:
  typedef V variable_t;
  typedef double value_type;
  typedef OptimisedPolynomial<variable_t> optimised_t;
  typedef CLNWrapper<precision> internal_value_t;
  typedef Monomial<variable_t, internal_value_t> monomial_t;
  typedef std::map<monomial_t, internal_value_t> coefficient_map_t;

private:
  util::LazyCopy<coefficient_map_t> coefficients;

  void addTerm(const internal_value_t& coefficient, const variable_t& variable, const std::size_t exponent)
  {
    (*coefficients)[monomial_t(variable, exponent)] += coefficient;
  }

  void addConstant(const internal_value_t& constant)
  {
    (*coefficients)[monomial_t()] += constant;
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

      if (currentIter->second == internal_value_t(0.0))
        coefficients->erase(currentIter);

      currentIter = nextIter;
    }
  }

  void addMonomial(const internal_value_t coefficient, const monomial_t& m)
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
  static const bool supports_abs = false;

  typedef detail::ValueMap<variable_t, internal_value_t> value_map;
  typedef typename coefficient_map_t::iterator iterator;
  typedef typename coefficient_map_t::const_iterator const_iterator;

  Polynomial()
  {
  }

  Polynomial(const Polynomial& p) : coefficients(p.coefficients)
  {
  }

  Polynomial(const value_type constant)
  {
    addConstant(constant);
    cleanZeros();
  }

  Polynomial(const variable_t& variable)
  {
    addTerm(1.0, variable, 1);
    cleanZeros();
  }

  Polynomial(const value_type coefficient, const variable_t& variable)
  {
    addTerm(coefficient, variable, 1);
    cleanZeros();
  }

  Polynomial(const variable_t& variable, const std::size_t exponent)
  {
    addTerm(1.0, variable, exponent);
    cleanZeros();
  }

  Polynomial(const value_type coefficient, const variable_t& variable, const std::size_t exponent)
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
    return *coefficients == *p.coefficients;
  }

  Polynomial& operator*=(const value_type x)
  {
    transformCoefficients(boost::bind(std::multiplies<internal_value_t>(), _1, x));
    cleanZeros();
    return *this;
  }

  Polynomial& operator/=(const value_type x)
  {
    transformCoefficients(boost::bind(std::divides<internal_value_t>(), _1, x));
    cleanZeros();
    return *this;
  }

  Polynomial& operator+=(const value_type x)
  {
    addConstant(x);
    cleanZeros();
    return *this;
  }

  Polynomial& operator-=(const value_type x)
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

  Polynomial pow(const int n) const
  {
    if (n<0)
      CFD_EXCEPTION("Polynomial cannot be raised to negative exponent.");

    Polynomial result(1);

    for(int i=0; i<n; ++i)
      result *= *this;

    return result;
  }

  std::size_t degree(const variable_t& variable) const
  {
    std::size_t result = 0;
    BOOST_FOREACH(const typename coefficient_map_t::value_type& cMapping, *coefficients)
    {
      result = std::max(result, cMapping.first.getExponent(variable));
    }
    return result;
  }

  Polynomial derivative(const variable_t& variable) const
  {
    Polynomial result;
    for(typename coefficient_map_t::const_iterator cIter(coefficients->begin()); cIter!=coefficients->end(); ++cIter)
    {
      const std::pair<internal_value_t, monomial_t> mDerivative(cIter->first.derivative(variable));
      result.addMonomial(cIter->second * mDerivative.first, mDerivative.second);
    }
  
    result.cleanZeros();
    return result;
  }

  Polynomial substituteValues(const value_map& valueMap) const
  {
    Polynomial result;
  
    for(typename coefficient_map_t::const_iterator cIter(coefficients->begin()); cIter!=coefficients->end(); ++cIter)
    {
      const std::pair<internal_value_t, monomial_t> mBound(cIter->first.substituteValues(valueMap));
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
  
      const bool renderCoefficient = cIter->second != internal_value_t(1.0) || cIter->first.isOne();
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

  void accept(NumericExpressionVisitor<variable_t>& v) const
  {
    BOOST_FOREACH(const typename coefficient_map_t::value_type& mapping, *coefficients)
    {
      const monomial_t& monomial = mapping.first;
      v.visitConstant(mapping.second);
      
      BOOST_FOREACH(const typename monomial_t::value_type& exponent, monomial)
      {
        v.visitVariable(exponent.first);
        v.visitExponent(exponent.second);
      }

      v.postProduct(1 + monomial.size());
    }

    v.postSummation(coefficients->size());
  }

  void swap(Polynomial& p)
  {
    coefficients.swap(p.coefficients);
  }
};

template<typename V>
Polynomial<V> pow(const Polynomial<V>& p, const int n)
{
  return p.pow(n);
}

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
