#ifndef SIMPLE_CFD_NUMERIC_EXCAFE_EXPRESSION_HPP
#define SIMPLE_CFD_NUMERIC_EXCAFE_EXPRESSION_HPP

#include <ostream>
#include <set>
#include <map>
#include <cassert>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include "expression.hpp"
#include "orthotope.hpp"
#include "expression_visitor.hpp"
#include "convert_expression.hpp"
#include "optimised_polynomial_fraction.hpp"
#include "excafe_mapper.hpp"
#include "excafe_value_map.hpp"
#include "symbol_collector.hpp"
#include <simple_cfd/exception.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/group.hpp>

namespace cfd
{

namespace 
{

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

   void visitConstant(const parent_t::integer_t& value)
   {
     visitor.visitConstant(value);
   }

   void visitConstant(const parent_t::float_t& value)
   {
     visitor.visitConstant(value);
   }

   void visitExponent(const int n)
   {
     visitor.visitExponent(n);
   }

   void visitAbsoluteValue()
   {
     visitor.visitAbsoluteValue();
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
  static const bool supports_abs = true;

  typedef double                                                value_type;
  typedef V                                                     variable_t;
  typedef ExcafeExpression<variable_t>                          optimised_t;
  typedef detail::ExcafeValueMap<variable_t, ExcafeExpression>  value_map;
  typedef numeric::Orthotope<variable_t, value_type>            region_t;

  friend class detail::ExcafeValueMap<variable_t, ExcafeExpression>;

private:
  typedef symbolic::Expr    expr_t;
  typedef symbolic::Symbol  symbol_t;
  
  expr_t expr;

  static symbol_t getSymbol(const variable_t& var)
  {
    detail::ExcafeMapper<variable_t>& mapper(detail::ExcafeMapper<variable_t>::instance());
    return mapper.getSymbol(var);
  }

  static variable_t getVariable(const symbol_t& s)
  {
    detail::ExcafeMapper<variable_t>& mapper(detail::ExcafeMapper<variable_t>::instance());
    return mapper.getKey(s);
  }

  explicit ExcafeExpression(const expr_t& e) : expr(e)
  {
  }

public:
  static ExcafeExpression group(const ExcafeExpression& e)
  {
    return ExcafeExpression(symbolic::Group(e.expr).clone());
  }

  ExcafeExpression()
  {
  }

  ExcafeExpression(const value_type s) : expr(s) 
  {
  }

  ExcafeExpression(const symbolic::Rational& s) : expr(s) 
  {
  }

  ExcafeExpression(const cln::cl_RA& s) : expr(s) 
  {
  }

  ExcafeExpression(const cln::cl_F& s) : expr(s) 
  {
  }

  ExcafeExpression(const variable_t& v) : expr(getSymbol(v)) 
  {
  }

  ExcafeExpression(const value_type c, const variable_t& v) : expr(expr_t(c) * getSymbol(v)) 
  {
  }

  ExcafeExpression(const variable_t& v, const std::size_t e) : expr(symbolic::pow(getSymbol(v), e)) 
  {
  }

  ExcafeExpression(const value_type c, const variable_t& v, const std::size_t e) : expr(c * symbolic::pow(getSymbol(v), e)) 
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

  ExcafeExpression integrate(const region_t& region) const
  {
    using namespace symbolic;

    Expr::region_t symbolicRegion;
    BOOST_FOREACH(const typename region_t::value_type& interval, region)
    {
      const Expr lowerBound = Float(interval.second.first).simplify();
      const Expr upperBound = Float(interval.second.second).simplify();

      if (is_a<Rational>(lowerBound) && is_a<Rational>(upperBound))
      {
        symbolicRegion.setInterval(getSymbol(interval.first), 
                                   convert_to<Rational>(lowerBound), 
                                   convert_to<Rational>(upperBound));
      }
      else
      {
        CFD_EXCEPTION("Can only integrate over rationally bounded regions.");
      }
    }

    const expr_t integrated = expr.integrate(symbolicRegion);
    return ExcafeExpression(integrated);
  }

  value_type operator()() const
  {
    value_map valueMap;
    return evaluate(valueMap);
  }

  value_type operator()(const value_type& v) const
  {
    const std::set<variable_t> variables = getVariables();
    assert(variables.size() <= 1);
    value_map valueMap;

    if (variables.size() == 1)
      valueMap.bind(*variables.begin(), v);

    return evaluate(valueMap);
  }

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
    return ExcafeExpression(expr.expand(true));
  }

  ExcafeExpression pow(const int n) const
  {
    return ExcafeExpression(symbolic::pow(expr, n));
  }
  
  ExcafeExpression abs() const
  {
    return ExcafeExpression(symbolic::abs(expr));
  }

  std::set<variable_t> getVariables() const
  {
    detail::SymbolCollector<variable_t> collector;
    accept(collector);
    return collector.getSymbols();
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
    const symbolic::Float evaluated = expr.eval(variableValues.getReference());
    return evaluated.toDouble();
  }

  boost::optional<symbolic::Rational> evaluateRational(const value_map& variableValues) const
  {
    using namespace symbolic;
    const expr_t evaluated = expr.subs(variableValues.getReference());

    if (is_exactly_a<Rational>(evaluated))
      return convert_to<Rational>(evaluated);
    else
      return boost::optional<Rational>();
  }
};

template<typename V>
std::ostream& operator<<(std::ostream& o, const ExcafeExpression<V>& e)
{
  e.write(o);
  return o;
}

template<typename V>
ExcafeExpression<V> pow(const ExcafeExpression<V>& e, const int n)
{
  return e.pow(n);
}

template<typename V>
ExcafeExpression<V> abs(const ExcafeExpression<V>& e)
{
  return e.abs();
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
