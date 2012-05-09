#ifndef EXCAFE_NUMERIC_OPTIMISED_POLYNOMIAL_HPP
#define EXCAFE_NUMERIC_OPTIMISED_POLYNOMIAL_HPP

#include <sstream>
#include <cmath>
#include <set>
#include <vector>
#include <utility>
#include <cstddef>
#include <cassert>
#include <ostream>
#include <boost/foreach.hpp>
#include "cast.hpp"
#include "numeric_fwd.hpp"
#include "polynomial.hpp"
#include "monomial.hpp"
#include "value_map.hpp"
#include <excafe/exception.hpp>
#include <excafe/util/hybrid_array.hpp>

namespace excafe
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
  static const std::size_t MAX_VALUES_STACK_BYTES = 512;
  static const std::size_t MAX_VALUES_STACK = MAX_VALUES_STACK_BYTES/sizeof(value_type);

  typedef std::pair<std::size_t, std::size_t> exponent_t;
  
  std::vector<exponent_t> exponents;
  std::vector<std::size_t> monomialLengths;
  std::vector<internal_value_t> monomialCoefficients;
  std::vector<variable_t> variables;

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

  value_type evaluate(const internal_value_t* const params) const
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

    return excafe::numeric_cast<value_type>(result);
  }

  template<typename element_t>
  static std::vector<element_t> asVector(const std::set<element_t>& s)
  {
    return std::vector<element_t>(s.begin(), s.end());
  }

public:
  OptimisedPolynomial()
  {
    monomialLengths.push_back(0);
    monomialCoefficients.push_back(0.0);
  }

  OptimisedPolynomial(const Polynomial<variable_t>& p) : variables(asVector(p.getVariables()))
  {
    p.checkConsistent();
  
    for(typename Polynomial<variable_t>::const_iterator mIter(p.begin()); mIter!=p.end(); ++mIter)
    {
      pushTerm(mIter->first, excafe::numeric_cast<internal_value_t>(mIter->second));
    }
  }

  std::set<variable_t> getVariables() const
  {
    return std::set<variable_t>(variables.begin(), variables.end());
  }

  value_type operator()() const
  {
    assert(variables.empty());
    return evaluate(NULL);
  }

  value_type operator()(const value_type a) const
  {
    /*
     HACK: We often find x cancels out in some of our univariate formulae, so we permit calling this
           where we have no variables.
    */
    assert(variables.size() < 2);
    internal_value_t paramData[1];
  
    if (variables.size() == 1)
      paramData[0] = a;

    return evaluate(paramData);
  }

  value_type operator()(const value_type a, const value_type b) const
  {
    assert(variables.size() == 2);
    internal_value_t paramData[2];
  
    paramData[0] = a;
    paramData[1] = b;
    return evaluate(paramData);
  }

  value_type operator()(const value_type a, const value_type b, const value_type c) const
  {
    assert(variables.size() == 3);
    internal_value_t paramData[3];
  
    paramData[0] = a;
    paramData[1] = b;
    paramData[2] = c;
    return evaluate(paramData);
  }

  bool isOne() const
  {
    return monomialLengths.size() == 1 && 
           monomialLengths[0] == 0 && 
           monomialCoefficients[0] == 1.0;
  }

  void write(std::ostream& out) const
  {
    assert(monomialCoefficients.size() == monomialLengths.size());

    std::size_t exponentIndex = 0;

    for(std::size_t monomialIndex=0; monomialIndex<monomialLengths.size(); ++monomialIndex)
    {
      const std::size_t newExponentIndex = exponentIndex + monomialLengths[monomialIndex];

      if (monomialIndex>0)
        out << " + ";

      out << monomialCoefficients[monomialIndex];

      if (exponentIndex != newExponentIndex)
        out << "*";

      for(; exponentIndex<newExponentIndex; ++exponentIndex)
      {
        out << variables[exponents[exponentIndex].first];
        const std::size_t exponent = exponents[exponentIndex].second;

        if (exponent > 1)
          out << "^" << exponent;
      }
    }

    if (monomialLengths.empty())
      out << 0.0;
  }

  value_type evaluate(const value_map& valueMap) const
  {
    typedef typename value_map::scalar_subst_map scalar_subst_map;

    util::HybridArray<internal_value_t, MAX_VALUES_STACK> paramData(variables.size());
    const scalar_subst_map& variableValues = valueMap.getScalarSubstitutions();
    typename scalar_subst_map::const_iterator varValIter = variableValues.begin();

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
        paramData[variableIndex] = excafe::numeric_cast<internal_value_t>(varValIter->second);
      }
    }

    return evaluate(paramData.get());
  }
};

template<typename V>
std::ostream& operator<<(std::ostream& o, const OptimisedPolynomial<V>& p)
{
  p.write(o);
  return o;
}

}

#endif
