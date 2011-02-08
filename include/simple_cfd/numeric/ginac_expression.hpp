#ifndef SIMPLE_CFD_NUMERIC_GINAC_EXPRESSION_HPP
#define SIMPLE_CFD_NUMERIC_GINAC_EXPRESSION_HPP

#include <ostream>
#include <sstream>
#include <set>
#include <map>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>
#include <ginac/ginac.h>
#include "ginac_mapper.hpp"
#include "ginac_value_map.hpp"
#include "expression.hpp"
#include "expression_visitor.hpp"
#include "convert_expression.hpp"
#include "optimised_polynomial_fraction.hpp"
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

template<typename V>
class GinacVisitorAdapter : public GiNaC::visitor, 
                            public GiNaC::mul::visitor, 
                            public GiNaC::add::visitor, 
                            public GiNaC::power::visitor, 
                            public GiNaC::symbol::visitor,
                            public GiNaC::numeric::visitor,
                            public GiNaC::basic::visitor
{
private:
  typedef V variable_t;
  NumericExpressionVisitor<variable_t>& visitor;

  void visitChildren(const GiNaC::expairseq& e)
  {
    for(std::size_t index=0; index<e.nops(); ++index)
      e.op(index).accept(*this);
  }

public:
  GinacVisitorAdapter(NumericExpressionVisitor<variable_t>& v) : visitor(v)
  {
  }

  void visit(const GiNaC::mul& b)
  {
    visitChildren(b);
    visitor.postProduct(b.nops());
  }

  void visit(const GiNaC::add& a)
  {
    visitChildren(a);
    visitor.postSummation(a.nops());
  }

  void visit(const GiNaC::symbol& s)
  {
    detail::GinacMapper<variable_t>& mapper(detail::GinacMapper<variable_t>::instance());
    visitor.visitVariable(mapper.getKey(s));
  }

  void visit(const GiNaC::numeric& n)
  {
    const cln::cl_R value = cln::realpart(n.to_cl_N());
    visitor.visitConstant(value);
  }

  void visit(const GiNaC::power& p)
  {
    p.op(0).accept(*this);
    const GiNaC::ex exponent = p.op(1);

    if (!GiNaC::is_a<GiNaC::numeric>(exponent))
    {
      CFD_EXCEPTION("NumericExpressionVisitor interface cannot handle non-numeric exponent.");
    }
    else
    {
      const GiNaC::numeric exponentNumeric = GiNaC::ex_to<GiNaC::numeric>(exponent);
      if (!exponentNumeric.is_integer())
      {
        CFD_EXCEPTION("NumericExpressionVisitor interface cannot handle non-integer exponent.");
      }
      else
      {
        visitor.visitExponent(exponentNumeric.to_long());
      }
    }
  }

  void visit(const GiNaC::basic& b)
  {
    std::ostringstream error;
    error << "Cannot handle " << b << " using NumericExpressionVisitor interface.";
    CFD_EXCEPTION(error.str());
  }
};

}

template<typename V>
class GinacExpression : public NumericExpression<V>,
                        boost::arithmetic<GinacExpression<V>, double,
                        boost::arithmetic<GinacExpression<V>
                        > >
{
public:
  typedef double                                              value_type;
  typedef V                                                   variable_t;
  typedef OptimisedPolynomialFraction<variable_t>             optimised_t;
  typedef detail::GinacValueMap<variable_t, GinacExpression>  value_map;
  friend class detail::GinacValueMap<variable_t, GinacExpression>;

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

public:
  GinacExpression() : expr(ginac_numeric_t(0.0))
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
    return *this;
  }

  GinacExpression& operator/=(const GinacExpression& e)
  {
    expr /= e.expr;
    return *this;
  }

  void accept(NumericExpressionVisitor<variable_t>& v) const
  {
    GinacVisitorAdapter<variable_t> adapter(v);
    expr.accept(adapter);
  }

  GinacExpression derivative(const variable_t& variable) const
  {
    return GinacExpression(expr.diff(getSymbol(variable)));
  }

  GinacExpression integrate(const variable_t& variable, const value_type& a, const value_type& b) const
  {
    return GinacExpression(GiNaC::integral(getSymbol(variable), a, b, expr).eval_integ().eval());
  }
  
  std::size_t degree(const variable_t& variable) const
  {
    return expr.degree(getSymbol(variable));
  }

  GinacExpression substituteValues(const value_map& valueMap) const
  {
    return GinacExpression(expr.subs(valueMap.getReference()));
  }

  optimised_t optimise() const
  {
    const GinacExpression normalised(expr.normal());
    const PolynomialFraction<variable_t> polyFraction = 
      cfd::detail::convert_expression< PolynomialFraction<variable_t> >(normalised);

    return optimised_t(polyFraction);
  }

  GinacExpression normalised() const
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
      const GiNaC::numeric numericValue = GiNaC::ex_to<GiNaC::numeric>(evaluated);
      return cfd::numeric_cast<value_type>(cln::cl_float(cln::realpart(numericValue.to_cl_N())));
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
