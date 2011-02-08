#ifndef SIMPLE_CFD_NUMERIC_EXCAFE_EXPRESSION_HPP
#define SIMPLE_CFD_NUMERIC_EXCAFE_EXPRESSION_HPP

#include <ostream>
#include <set>
#include <map>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>
#include "expression.hpp"
#include "expression_visitor.hpp"
#include "convert_expression.hpp"
#include "optimised_polynomial_fraction.hpp"
#include "excafe_mapper.hpp"
#include "excafe_value_map.hpp"
#include <simple_cfd/exception.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/symbol.hpp>

namespace cfd
{

namespace 
{

class ExcafeSymbolCollector : public symbolic::Visitor, 
                              public symbolic::Symbol::Visitor
{
private:
  std::set<symbolic::Symbol> symbols;

public:
  void visit(const symbolic::Symbol& s) 
  {
    symbols.insert(s);
  }

  std::set<symbolic::Symbol> getSymbols() const
  {
    return symbols;
  }
};

template<typename V>
class ExcafeVisitorAdapter : public NumericExpressionVisitor<symbolic::Symbol>
{
private:
   typedef NumericExpressionVisitor<symbolic::Symbol> parent_t;
   typedef V variable_t;
   NumericExpressionVisitor<variable_t>& visitor;

public:
   ExcafeVisitorAdapter(NumericExpressionVisitor<variable_t>& v) : visitor(v)
   {
   }

   void visitConstant(const parent_t::value_t& value)
   {
     visitor.visitConstant(value);
   }

   void visitExponent(const int n)
   {
     visitor.visitExponent(n);
   }

   void postSummation(const std::size_t n)
   {
     visitor.postSummation(n);
   }

   void postProduct(const std::size_t n)
   {
     visitor.postProduct(n);
   }

   void visitVariable(const symbolic::Symbol& v)
   {
     detail::ExcafeMapper<variable_t>& mapper(detail::ExcafeMapper<variable_t>::instance());
     visitor.visitVariable(mapper.getKey(v));
   }
};

}
template<typename V>
class ExcafeExpression : public NumericExpression<V>,
                        boost::arithmetic<ExcafeExpression<V>, double,
                        boost::arithmetic<ExcafeExpression<V>
                        > >
{
public:
  typedef double                                                value_type;
  typedef V                                                     variable_t;
  typedef ExcafeExpression<variable_t>                          optimised_t;
  typedef detail::ExcafeValueMap<variable_t, ExcafeExpression>  value_map;
  friend class detail::ExcafeValueMap<variable_t, ExcafeExpression>;

private:
  typedef symbolic::Expr    expr_t;
  typedef symbolic::Symbol  symbol_t;
  typedef symbolic::Float   numeric_t;
  
  expr_t expr;

  static symbol_t getSymbol(const variable_t& var)
  {
    detail::ExcafeMapper<variable_t>& mapper(detail::ExcafeMapper<variable_t>::instance());
    return mapper.getExcafeSymbol(var);
  }

  static variable_t getVariable(const symbol_t& s)
  {
    detail::ExcafeMapper<variable_t>& mapper(detail::ExcafeMapper<variable_t>::instance());
    return mapper.getKey(s);
  }

  ExcafeExpression(const expr_t& e) : expr(e)
  {
  }

public:
  ExcafeExpression() : expr(numeric_t(0.0))
  {
  }

  ExcafeExpression(const value_type s) : expr(numeric_t(s)) 
  {
  }

  ExcafeExpression(const variable_t& v) : expr(getSymbol(v)) 
  {
  }

  ExcafeExpression(const value_type c, const variable_t& v) : expr(expr_t(c) * getSymbol(v)) 
  {
  }

  ExcafeExpression(const variable_t& v, const std::size_t e) : expr(pow(getSymbol(v), e)) 
  {
  }

  ExcafeExpression(const value_type c, const variable_t& v, const std::size_t e) : expr(c * pow(getSymbol(v), e)) 
  {
  }

  ExcafeExpression& operator=(const ExcafeExpression& e)
  {
    expr = e.expr;
    return *this;
  }

  ExcafeExpression& operator+=(const value_type s)
  {
    expr += s;
    return *this;
  }

  ExcafeExpression& operator-=(const value_type s)
  {
    expr -= s;
    return *this;
  }

  ExcafeExpression& operator*=(const value_type s)
  {
    expr *= s;
    return *this;
  }

  ExcafeExpression& operator/=(const value_type s)
  {
    expr /= s;
    return *this;
  }

  ExcafeExpression operator-() const
  {
    return ExcafeExpression(-expr);
  }

  ExcafeExpression& operator+=(const ExcafeExpression& e)
  {
    expr += e.expr;
    return *this;
  }

  ExcafeExpression& operator-=(const ExcafeExpression& e)
  {
    expr -= e.expr;
    return *this;
  }

  ExcafeExpression& operator*=(const ExcafeExpression& e)
  {
    expr *= e.expr;
    return *this;
  }

  ExcafeExpression& operator/=(const ExcafeExpression& e)
  {
    expr /= e.expr;
    return *this;
  }

  void accept(NumericExpressionVisitor<variable_t>& v) const
  {
    ExcafeVisitorAdapter<variable_t> adapter(v);
    expr.accept(adapter);
  }

  ExcafeExpression derivative(const variable_t& variable) const
  {
    return ExcafeExpression(expr.derivative(getSymbol(variable)));
  }

  ExcafeExpression integrate(const variable_t& variable, const value_type a, const value_type b) const
  {
    return expr.integrate(getSymbol(variable), a, b);
  }

/*
  std::size_t degree(const variable_t& variable) const
  {
    return expr.degree(getSymbol(variable));
  }
*/

  ExcafeExpression substituteValues(const value_map& valueMap) const
  {
    return ExcafeExpression(expr.subs(valueMap.getReference()));
  }

  optimised_t optimise() const
  {
    return ExcafeExpression(expr.simplify());
  }

  ExcafeExpression normalised() const
  {
    return *this;
  }

  std::set<variable_t> getVariables() const
  {
    ExcafeSymbolCollector collector;
    expr.traverse(collector);

    const std::set<symbolic::Symbol> symbols(collector.getSymbols());
    std::set<variable_t> variables;

    BOOST_FOREACH(const symbolic::Symbol& s, symbols)
    {
      variables.insert(getVariable(s));
    }

    return variables;
  }

  void write(std::ostream& o) const
  {
    o << expr;
  }

  void swap(ExcafeExpression& e)
  {
    expr.swap(e.expr);
  }

  value_type evaluate(const value_map& variableValues) const
  {
    const expr_t evaluated = expr.subs(variableValues.getReference());
    const symbolic::Float* const f = dynamic_cast<const symbolic::Float*>(&expr.internal());

    if (f!=NULL)
    {
      return f->toDouble();
    }
    else
    {
      CFD_EXCEPTION("Evaluation of Excafe expression failed to produce a numeric value.");
    }
  }
};

template<typename V>
std::ostream& operator<<(std::ostream& o, const ExcafeExpression<V>& e)
{
  e.write(o);
  return o;
}

}

namespace std
{

template<typename V>
void swap(cfd::ExcafeExpression<V>& a, cfd::ExcafeExpression<V>& b)
{
  a.swap(b);
}

}

#endif
