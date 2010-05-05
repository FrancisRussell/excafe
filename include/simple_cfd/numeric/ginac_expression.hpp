#ifndef SIMPLE_CFD_NUMERIC_GINAC_EXPRESSION_HPP
#define SIMPLE_CFD_NUMERIC_GINAC_EXPRESSION_HPP

#include <ostream>
#include <boost/operators.hpp>
#include <ginac/ginac.h>
#include "ginac_mapper.hpp"

namespace cfd
{

template<typename V>
class GinacExpression : boost::addable<GinacExpression<V>, double,
                        boost::subtractable<GinacExpression<V>, double,
                        boost::multipliable<GinacExpression<V>, double,
                        boost::dividable<GinacExpression<V>, double,
                        boost::addable<GinacExpression<V>,
                        boost::subtractable<GinacExpression<V>,
                        boost::dividable<GinacExpression<V>,
                        boost::multipliable<GinacExpression<V>
                        > > > > > > > >
{
public:
  typedef V variable_t;

private:
  typedef GiNaC::ex      ginac_expr_t;
  typedef GiNaC::symbol  ginac_symbol_t;
  typedef GiNaC::numeric ginac_numeric_t;
  
  ginac_expr_t expr;

  ginac_symbol_t getSymbol(const variable_t& var) const
  {
    detail::GinacMapper<variable_t>& mapper(detail::GinacMapper<variable_t>::instance());
    return mapper.getSymbol(var);
  }

  GinacExpression(const ginac_expr_t& e) : expr(e)
  {
  }

public:
  GinacExpression()
  {
  }

  GinacExpression(const double s) : expr(ginac_numeric_t(s)) 
  {
  }

  GinacExpression(const variable_t& v) : expr(getSymbol(v)) 
  {
  }

  GinacExpression(const double c, const variable_t& v) : expr(c * getSymbol(v)) 
  {
  }

  GinacExpression(const variable_t& v, const std::size_t e) : expr(pow(getSymbol(v), e)) 
  {
  }

  GinacExpression(const double c, const variable_t& v, const std::size_t e) : expr(c * pow(getSymbol(v), e)) 
  {
  }

  GinacExpression& operator=(const GinacExpression& e)
  {
    expr = e.expr;
    return *this;
  }

  GinacExpression& operator+=(const double s)
  {
    expr += s;
    return *this;
  }

  GinacExpression& operator-=(const double s)
  {
    expr -= s;
    return *this;
  }

  GinacExpression& operator*=(const double s)
  {
    expr *= s;
    return *this;
  }

  GinacExpression& operator/=(const double s)
  {
    expr /= s;
    return *this;
  }

  GinacExpression operator-() const
  {
    return GinacExpression(-expr);
  }

  GinacExpression& operator+=(const GinacExpression& e)
  {
    expr += e.expr;
    return *this;
  }

  GinacExpression& operator-=(const GinacExpression& e)
  {
    expr -= e.expr;
    return *this;
  }

  GinacExpression& operator*=(const GinacExpression& e)
  {
    expr *= e.expr;
    return *this;
  }

  GinacExpression& operator/=(const GinacExpression& e)
  {
    expr /= e.expr;
    return *this;
  }

  GinacExpression derivative(const variable_t& variable) const
  {
    return GinacExpression(expr.diff(getSymbol(variable)));
  }

  void write(std::ostream& o) const
  {
    o << GiNaC::dflt << expr;
  }

  void swap(GinacExpression& e)
  {
    expr.swap(e.expr);
  }
};

template<typename V>
std::ostream& operator<<(std::ostream& o, const GinacExpression<V>& e)
{
  e.write(o);
  return o;
}

}

namespace std
{

template<typename V>
void swap(cfd::GinacExpression<V>& a, cfd::GinacExpression<V>& b)
{
  a.swap(b);
}

}

#endif
