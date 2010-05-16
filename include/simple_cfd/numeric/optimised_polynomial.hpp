#ifndef SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_HPP

#include <sstream>
#include <numeric>
#include <cmath>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <cstddef>
#include <cassert>
#include <boost/foreach.hpp>
#include "cast.hpp"
#include "numeric_fwd.hpp"
#include "polynomial.hpp"
#include "monomial.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace 
{

template<typename T>
struct Pow
{
  typedef T value_type;

  inline value_type operator()(const value_type value, const std::size_t exponent) const
  {
    return std::pow(value, exponent);
  }
};

}

template<typename V>
class OptimisedPolynomial
{
public:
  typedef V variable_t;
  typedef double value_type;
  typedef typename Polynomial<variable_t>::value_map value_map;

private:
  typedef std::vector< std::pair<std::vector<std::size_t>, value_type> > coefficient_vec_t;
  std::set<variable_t> variables;
  coefficient_vec_t coefficients;

  // A slightly hacky solution to avoid dynamic memory allocation when evaluating.
  mutable std::vector<value_type> paramData;

  template<typename monomial_t>
  std::vector<std::size_t> buildExponentVector(const monomial_t& m) const
  {
    std::vector<std::size_t> exponents;
    exponents.reserve(variables.size());
  
    BOOST_FOREACH(const variable_t& var, variables)
    {
      exponents.push_back(m.getExponent(var));
    }
  
    return exponents;
  }

  value_type evaluate(const std::vector<value_type>& params) const
  {
    value_type result = 0.0;
  
    for(coefficient_vec_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
    {
      result += std::inner_product(params.begin(), params.end(), cIter->first.begin(), cIter->second,
        std::multiplies<value_type>(), Pow<value_type>());
    }
  
    return result;
  }

public:
  OptimisedPolynomial()
  {
  }

  OptimisedPolynomial(const Polynomial<variable_t>& p) : variables(p.getVariables()),
    paramData(variables.size())
  {
    p.checkConsistent();
  
    for(typename Polynomial<variable_t>::const_iterator mIter(p.begin()); mIter!=p.end(); ++mIter)
      coefficients.push_back(std::make_pair(buildExponentVector(mIter->first),
        cfd::numeric_cast<value_type>(mIter->second)));
  }

  std::set<variable_t> getVariables() const
  {
    return variables;
  }

  value_type operator()() const
  {
    assert(variables.empty());
    return evaluate(paramData);
  }

  value_type operator()(const value_type a) const
  {
    assert(variables.size() == 1);
  
    paramData[0] = a;
    return evaluate(paramData);
  }

  value_type operator()(const value_type a, const value_type b) const
  {
    assert(variables.size() == 2);
  
    paramData[0] = a;
    paramData[1] = b;
    return evaluate(paramData);
  }

  value_type operator()(const value_type a, const value_type b, const value_type c) const
  {
    assert(variables.size() == 3);
  
    paramData[0] = a;
    paramData[1] = b;
    paramData[2] = c;
    return evaluate(paramData);
  }

  value_type evaluate(const value_map& valueMap) const
  {
    typedef typename value_map::internal_map_t internal_value_map_t;
    const internal_value_map_t& variableValues = valueMap.getReference();

    std::size_t variableIndex = 0;
    typename internal_value_map_t::const_iterator varValIter = variableValues.begin();

    BOOST_FOREACH(const variable_t& v, variables)
    {
      while (varValIter != variableValues.end() && varValIter->first != v)
        ++varValIter;

      if (varValIter == variableValues.end())
      {
        std::ostringstream error;
        error << "Missing variable binding when evaluating OptimisedPolynomial: " << v << ".";
        CFD_EXCEPTION(error.str());
      }
      else
      {
        paramData[variableIndex] = cfd::numeric_cast<value_type>(varValIter->second);
      }
      ++variableIndex;
    }

    return evaluate(paramData);
  }
};

}

#endif
