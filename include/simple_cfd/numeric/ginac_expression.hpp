#ifndef SIMPLE_CFD_NUMERIC_GINAC_EXPRESSION_HPP
#define SIMPLE_CFD_NUMERIC_GINAC_EXPRESSION_HPP

#include <ostream>
#include <set>
#include <map>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>
#include <ginac/ginac.h>
#include "ginac_mapper.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace 
{

class GinacSymbolCollector : public GiNaC::visitor, public GiNaC::symbol::visitor
{
private:
  std::set<GiNaC::symbol> symbols;

public:
  void visit(const GiNaC::symbol& s) 
  {
    symbols.insert(s);
  }

  std::set<GiNaC::symbol> getSymbols() const
  {
    return symbols;
  }
};

}

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
  typedef double value_type;
  typedef V variable_t;
  typedef GinacExpression<variable_t> optimised_t;
  typedef std::map<variable_t, value_type>   value_map_t;

private:
  typedef GiNaC::ex      ginac_expr_t;
  typedef GiNaC::symbol  ginac_symbol_t;
  typedef GiNaC::numeric ginac_numeric_t;
  
  ginac_expr_t expr;

  ginac_symbol_t getSymbol(const variable_t& var) const
  {
    detail::GinacMapper<variable_t>& mapper(detail::GinacMapper<variable_t>::instance());
    return mapper.getGiNaCSymbol(var);
  }

  variable_t getVariable(const ginac_symbol_t& s) const
  {
    detail::GinacMapper<variable_t>& mapper(detail::GinacMapper<variable_t>::instance());
    return mapper.getKey(s);
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

  GinacExpression substituteValue(const variable_t& variable, const double value) const
  {
    return GinacExpression(expr.subs(getSymbol(variable) == value));
  }

  optimised_t optimise() const
  {
    return GinacExpression(expr.normal());
  }

  std::set<variable_t> getVariables() const
  {
    GinacSymbolCollector collector;
    expr.accept(collector);

    const std::set<GiNaC::symbol> ginacSymbols(collector.getSymbols());
    std::set<variable_t> variables;

    BOOST_FOREACH(const GiNaC::symbol& s, ginacSymbols)
    {
      variables.insert(getVariable(s));
    }

    return variables;
  }

  void write(std::ostream& o) const
  {
    o << GiNaC::dflt << expr;
  }

  void swap(GinacExpression& e)
  {
    expr.swap(e.expr);
  }

  value_type evaluate(const value_map_t& variableValues) const
  {
    typedef GiNaC::exmap ginac_exmap_t;
    ginac_exmap_t ginacMap;

    BOOST_FOREACH(const typename value_map_t::value_type& mapping, variableValues)
    {
      ginacMap.insert(ginac_exmap_t::value_type(getSymbol(mapping.first), ginac_numeric_t(mapping.second)));
    }

    const ginac_expr_t evaluated = expr.subs(ginacMap);

    if (GiNaC::is_a<GiNaC::numeric>(evaluated))
    {
      return GiNaC::ex_to<GiNaC::numeric>(evaluated).to_double();
    }
    else
    {
      CFD_EXCEPTION("Evaluation of GiNaC expression failed to produce a numeric value.");
    }
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
