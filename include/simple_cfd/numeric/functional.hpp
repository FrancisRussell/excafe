#ifndef SIMPLE_CFD_NUMERIC_FUNCTIONAL_HPP
#define SIMPLE_CFD_NUMERIC_FUNCTIONAL_HPP

#include <set>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include "numeric_fwd.hpp"

namespace cfd
{

template<typename T>
class PolynomialDifferentiator
{
private:
  typedef T value_t;
  value_t dx;

public:
  PolynomialDifferentiator(const value_t& _dx) : dx(_dx)
  {
  }

  template<typename P>
  P operator()(const P& p) const
  {
    BOOST_STATIC_ASSERT((boost::is_same<P, Polynomial<value_t> >::value || 
      boost::is_same<P, PolynomialFraction<value_t> >::value));
    return p.derivative(dx);
  }
};

template<typename T>
class PolynomialOptimiser
{
public:
  typedef T value_type;
  typedef OptimisedPolynomial<value_type> result_type;

  result_type operator()(const Polynomial<value_type>& p) const
  {
    return p;
  }
};

template<typename T>
class PolynomialFractionOptimiser
{
public:
  typedef T value_type;
  typedef OptimisedPolynomialFraction<value_type> result_type;

  result_type operator()(const PolynomialFraction<value_type>& p) const
  {
    return p;
  }
};

template<typename T>
class PolynomialVariableCollector
{
private:
  typedef T value_type;
  std::set<value_type> variables;

public:
  template<typename Poly>
  void operator()(const Poly& p)
  {
    const std::set<value_type> vars(p.getVariables());
    variables.insert(vars.begin(), vars.end());
  }

  std::set<value_type> getVariables() const
  {
    return variables;
  }
};


}

#endif
