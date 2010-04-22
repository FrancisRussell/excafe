#ifndef SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_HPP

#include "numeric_fwd.hpp"
#include <numeric/polynomial.hpp>
#include <numeric/monomial.hpp>
#include <numeric>
#include <set>
#include <vector>
#include <utility>
#include <cstddef>
#include <cassert>

namespace cfd
{

namespace 
{

struct Pow
{
  inline double operator()(const double value, const std::size_t exponent) const
  {
    double result=1.0;

    for(std::size_t power=0; power<exponent; ++power)
      result *= value;

    return result;
  }
};

}

template<typename V>
class OptimisedPolynomial
{
public:
  typedef V variable_t;

private:
  typedef std::vector< std::pair<std::vector<std::size_t>, double> > coefficient_vec_t;
  std::set<variable_t> variables;
  coefficient_vec_t coefficients;
  mutable std::vector<double> paramData;

  std::vector<std::size_t> buildExponentVector(const Monomial<variable_t>& m) const
  {
    std::vector<std::size_t> exponents;
  
    for(typename std::set<variable_t>::const_iterator varIter(variables.begin()); varIter!=variables.end(); ++varIter)
      exponents.push_back(m.getExponent(*varIter));
  
    return exponents;
  }

  double evaluate(const std::vector<double>& params) const
  {
    double result = 0.0;
  
    for(coefficient_vec_t::const_iterator cIter(coefficients.begin()); cIter!=coefficients.end(); ++cIter)
      result += std::inner_product(params.begin(), params.end(), cIter->first.begin(), cIter->second, std::multiplies<double>(), Pow());
  
    return result;
  }

public:
  OptimisedPolynomial()
  {
  }

  OptimisedPolynomial(const Polynomial<variable_t>& p) : variables(p.getIndependentVariables()),
  paramData(variables.size())
  {
    p.checkConsistent();
  
    for(typename Polynomial<variable_t>::const_iterator mIter(p.begin()); mIter!=p.end(); ++mIter)
      coefficients.push_back(std::make_pair(buildExponentVector(mIter->first), mIter->second));
  }

  double operator()() const
  {
    assert(variables.empty());
    return evaluate(paramData);
  }

  double operator()(const double a) const
  {
    assert(variables.size() == 1);
  
    paramData[0] = a;
    return evaluate(paramData);
  }

  double operator()(const double a, const double b) const
  {
    assert(variables.size() == 2);
  
    paramData[0] = a;
    paramData[1] = b;
    return evaluate(paramData);
  }

  double operator()(const double a, const double b, const double c) const
  {
    assert(variables.size() == 3);
  
    paramData[0] = a;
    paramData[1] = b;
    paramData[2] = c;
    return evaluate(paramData);
  }

};

}

#endif
