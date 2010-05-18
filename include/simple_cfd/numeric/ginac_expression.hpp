#ifndef SIMPLE_CFD_NUMERIC_GINAC_EXPRESSION_HPP
#define SIMPLE_CFD_NUMERIC_GINAC_EXPRESSION_HPP

#include <ostream>
#include <set>
#include <map>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>
#include <ginac/ginac.h>
#include "ginac_mapper.hpp"
#include "ginac_value_map.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace 
{

class GinacSymbolCollector : public GiNaC::visitor, public GiNaC::symbol::visitor
{
private:
  std::set<GiNaC::symbol, GiNaC::ex_is_less> symbols;

public:
  void visit(const GiNaC::symbol& s) 
  {
    symbols.insert(s);
  }

  std::set<GiNaC::symbol, GiNaC::ex_is_less> getSymbols() const
  {
    return symbols;
  }
};

}

template<typename V>
class GinacExpression : boost::arithmetic<GinacExpression<V>, double,
                        boost::arithmetic<GinacExpression<V>
                        > >
{
public:
  typedef double                             value_type;
  typedef V                                  variable_t;
  typedef GinacExpression<variable_t>        optimised_t;
  typedef detail::GinacValueMap<variable_t>  value_map;

private:
  typedef GiNaC::ex      ginac_expr_t;
  typedef GiNaC::symbol  ginac_symbol_t;
  typedef GiNaC::numeric ginac_numeric_t;
  
  ginac_expr_t expr;

  static ginac_symbol_t getSymbol(const variable_t& var)
  {
    detail::GinacMapper<variable_t>& mapper(detail::GinacMapper<variable_t>::instance());
    return mapper.getGiNaCSymbol(var);
  }

  static variable_t getVariable(const ginac_symbol_t& s)
  {
    detail::GinacMapper<variable_t>& mapper(detail::GinacMapper<variable_t>::instance());
    return mapper.getKey(s);
  }

  GinacExpression(const ginac_expr_t& e) : expr(e)
  {
  }

  void simplify()
  {
    expr = expr.normal();
  }

public:
  GinacExpression()
  {
  }

  GinacExpression(const value_type s) : expr(ginac_numeric_t(s)) 
  {
  }

  GinacExpression(const variable_t& v) : expr(getSymbol(v)) 
  {
  }

  GinacExpression(const value_type c, const variable_t& v) : expr(c * getSymbol(v)) 
  {
  }

  GinacExpression(const variable_t& v, const std::size_t e) : expr(pow(getSymbol(v), e)) 
  {
  }

  GinacExpression(const value_type c, const variable_t& v, const std::size_t e) : expr(c * pow(getSymbol(v), e)) 
  {
  }

  GinacExpression& operator=(const GinacExpression& e)
  {
    expr = e.expr;
    return *this;
  }

  GinacExpression& operator+=(const value_type s)
  {
    expr += s;
    return *this;
  }

  GinacExpression& operator-=(const value_type s)
  {
    expr -= s;
    return *this;
  }

  GinacExpression& operator*=(const value_type s)
  {
    expr *= s;
    return *this;
  }

  GinacExpression& operator/=(const value_type s)
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
    simplify();
    return *this;
  }

  GinacExpression& operator/=(const GinacExpression& e)
  {
    expr /= e.expr;
    simplify();
    return *this;
  }

  GinacExpression derivative(const variable_t& variable) const
  {
    return GinacExpression(expr.diff(getSymbol(variable)));
  }

  GinacExpression substituteValues(const value_map& valueMap) const
  {
    return GinacExpression(expr.subs(valueMap.getReference()));
  }

  optimised_t optimise() const
  {
    return GinacExpression(expr.normal());
  }

  std::set<variable_t> getVariables() const
  {
    GinacSymbolCollector collector;
    expr.traverse(collector);

    const std::set<GiNaC::symbol, GiNaC::ex_is_less> ginacSymbols(collector.getSymbols());
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

  value_type evaluate(const value_map& variableValues) const
  {
    const ginac_expr_t evaluated = GiNaC::evalf(expr.subs(variableValues.getReference()));

    if (GiNaC::is_a<GiNaC::numeric>(evaluated))
    {
      return cfd::numeric_cast<value_type>(GiNaC::ex_to<GiNaC::numeric>(evaluated).to_cl_N());
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
