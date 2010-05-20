#ifndef SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_HPP

#include <sstream>
#include <cmath>
#include <set>
#include <vector>
#include <utility>
#include <cstddef>
#include <cassert>
#include <boost/foreach.hpp>
#include "cast.hpp"
#include "numeric_fwd.hpp"
#include "polynomial.hpp"
#include "monomial.hpp"
#include "value_map.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

template<typename V>
class OptimisedPolynomial
{
private:
  typedef double internal_value_t;

public:
  typedef V variable_t;
  typedef double value_type;
  typedef detail::ValueMap<variable_t, internal_value_t> value_map;

private:
  typedef std::pair<std::size_t, std::size_t> exponent_t;
  
  std::vector<exponent_t> exponents;
  std::vector<std::size_t> monomialLengths;
  std::vector<internal_value_t> monomialCoefficients;
  std::vector<variable_t> variables;

  // A slightly hacky solution to avoid dynamic memory allocation when evaluating.
  mutable std::vector<internal_value_t> paramData;

  template<typename M>
  void pushTerm(const M& monomial, const internal_value_t& coefficient)
  {
    pushMonomial(monomial);
    monomialCoefficients.push_back(coefficient);
  }

  template<typename M>
  void pushMonomial(const M& m)
  {
    std::size_t monomialLength = 0;

    for(std::size_t varIndex=0; varIndex < variables.size(); ++varIndex)
    {
      const variable_t& var = variables[varIndex];
      const std::size_t exponent = m.getExponent(var);
      if (exponent > 0)
      {
        exponents.push_back(std::make_pair(varIndex, exponent));
        ++monomialLength;
      }
    }
  
    monomialLengths.push_back(monomialLength);
  }

  value_type evaluate(const std::vector<internal_value_t>& params) const
  {
    assert(monomialCoefficients.size() == monomialLengths.size());

    internal_value_t result = 0.0;
    std::size_t exponentIndex = 0;

    for(std::size_t monomialIndex=0; monomialIndex<monomialLengths.size(); ++monomialIndex)
    {
      const std::size_t newExponentIndex = exponentIndex + monomialLengths[monomialIndex];
      internal_value_t monomialValue = 1.0;

      for(; exponentIndex<newExponentIndex; ++exponentIndex)
      {
        monomialValue *= pow(params[exponents[exponentIndex].first], exponents[exponentIndex].second);
      }

      result += monomialValue * monomialCoefficients[monomialIndex];
    }

    return cfd::numeric_cast<value_type>(result);
  }

  template<typename element_t>
  static std::vector<element_t> asVector(const std::set<element_t>& s)
  {
    return std::vector<element_t>(s.begin(), s.end());
  }

public:
  OptimisedPolynomial()
  {
  }

  OptimisedPolynomial(const Polynomial<variable_t>& p) : variables(asVector(p.getVariables())),
    paramData(variables.size())
  {
    p.checkConsistent();
  
    for(typename Polynomial<variable_t>::const_iterator mIter(p.begin()); mIter!=p.end(); ++mIter)
    {
      pushTerm(mIter->first, cfd::numeric_cast<internal_value_t>(mIter->second));
    }
  }

  std::set<variable_t> getVariables() const
  {
    return std::set<variable_t>(variables.begin(), variables.end());
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
    typename internal_value_map_t::const_iterator varValIter = variableValues.begin();

    for(std::size_t variableIndex=0; variableIndex < variables.size(); ++variableIndex)
    {
      const variable_t& v = variables[variableIndex];

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
        paramData[variableIndex] = cfd::numeric_cast<internal_value_t>(varValIter->second);
      }
    }

    return evaluate(paramData);
  }
};

}

#endif
