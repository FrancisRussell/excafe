#ifndef SIMPLE_CFD_NUMERIC_FUNCTIONAL_HPP
#define SIMPLE_CFD_NUMERIC_FUNCTIONAL_HPP

#include <set>
#include <map>
#include "numeric_fwd.hpp"

namespace cfd
{

template<typename P>
class PolynomialDifferentiator
{
private:
  typedef P polynomial_t;
  typedef typename polynomial_t::variable_t variable_t;
  variable_t dx;

public:
  PolynomialDifferentiator(const variable_t& _dx) : dx(_dx)
  {
  }

  polynomial_t operator()(const polynomial_t& p) const
  {
    return p.derivative(dx);
  }
};

template<typename P>
class PolynomialOptimiser
{
public:
  typedef P polynomial_t;
  typedef typename polynomial_t::optimised_t result_type;

  result_type operator()(const polynomial_t& p) const
  {
    return p.optimise();
  }
};

template<typename P>
class PolynomialVariableCollector
{
private:
  typedef P polynomial_t;
  typedef typename polynomial_t::variable_t variable_t;
  std::set<variable_t> variables;

public:
  void operator()(const polynomial_t& p)
  {
    const std::set<variable_t> vars(p.getVariables());
    variables.insert(vars.begin(), vars.end());
  }

  std::set<variable_t> getVariables() const
  {
    return variables;
  }
};

template<typename P>
class PolynomialEvaluator
{
private:
  typedef P polynomial_t;
  typedef typename polynomial_t::variable_t  variable_t;
  typedef typename polynomial_t::value_type  value_type;
  typedef typename polynomial_t::value_map   value_map;

  value_map variableValues;

public:
  typedef value_type result_type;

  PolynomialEvaluator(const value_map& _variableValues) :
    variableValues(_variableValues)
  {
  }

  result_type operator()(const polynomial_t& p) const
  {
    return p.evaluate(variableValues);
  }
};

}

#endif
