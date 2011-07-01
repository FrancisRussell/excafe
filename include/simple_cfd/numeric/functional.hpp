#ifndef SIMPLE_CFD_NUMERIC_FUNCTIONAL_HPP
#define SIMPLE_CFD_NUMERIC_FUNCTIONAL_HPP

#include <set>
#include <map>
#include "numeric_fwd.hpp"

namespace cfd
{

template<typename E>
class ExpressionDifferentiator
{
private:
  typedef E expression_t;
  typedef typename expression_t::variable_t variable_t;
  variable_t x;

public:
  ExpressionDifferentiator(const variable_t& _x) : x(_x)
  {
  }

  expression_t operator()(const expression_t& e) const
  {
    return e.derivative(x);
  }
};

template<typename E>
class ExpressionOptimiser
{
public:
  typedef E expression_t;
  typedef typename expression_t::optimised_t result_type;

  result_type operator()(const expression_t& e) const
  {
    return e.optimise();
  }
};

template<typename E>
class ExpressionVariableCollector
{
private:
  typedef E expression_t;
  typedef typename expression_t::variable_t variable_t;
  std::set<variable_t> variables;

public:
  void operator()(const expression_t& e)
  {
    const std::set<variable_t> vars(e.getVariables());
    variables.insert(vars.begin(), vars.end());
  }

  std::set<variable_t> getVariables() const
  {
    return variables;
  }
};

template<typename E>
class ExpressionEvaluator
{
private:
  typedef E expression_t;
  typedef typename expression_t::variable_t  variable_t;
  typedef typename expression_t::value_type  value_type;
  typedef typename expression_t::value_map   value_map;

  value_map variableValues;

public:
  typedef value_type result_type;

  ExpressionEvaluator(const value_map& _variableValues) :
    variableValues(_variableValues)
  {
  }

  result_type operator()(const expression_t& e) const
  {
    return e.evaluate(variableValues);
  }
};

}

#endif
