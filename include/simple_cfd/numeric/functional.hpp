#ifndef SIMPLE_CFD_NUMERIC_FUNCTIONAL_HPP
#define SIMPLE_CFD_NUMERIC_FUNCTIONAL_HPP

#include "numeric_fwd.hpp"
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

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

}

#endif
