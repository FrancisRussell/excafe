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
#include <simple_cfd/symbolic/number.hpp>
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
/*
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
    const cln::cl_F value = cln::cl_float(cln::realpart(n.to_cl_N()));
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
*/

}
template<typename V>
class ExcafeExpression /*: public NumericExpression<V>,
                        boost::arithmetic<ExcafeExpression<V>, double,
                        boost::arithmetic<ExcafeExpression<V>
                        > >*/
{
public:
  typedef double                                                value_type;
  typedef V                                                     variable_t;
//  typedef OptimisedPolynomialFraction<variable_t>             optimised_t;
  typedef detail::ExcafeValueMap<variable_t, ExcafeExpression>  value_map;
//  friend class detail::GinacValueMap<variable_t, ExcafeExpression>;

private:
  typedef symbolic::Expr    expr_t;
  typedef symbolic::Symbol  symbol_t;
  typedef symbolic::Number  numeric_t;
  
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

  ExcafeExpression(const value_type c, const variable_t& v) : expr(c * getSymbol(v)) 
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
/*
  void accept(NumericExpressionVisitor<variable_t>& v) const
  {
    GinacVisitorAdapter<variable_t> adapter(v);
    expr.accept(adapter);
  }
*/
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

/*
  optimised_t optimise() const
  {
    const ExcafeExpression normalised(expr.normal());
    const PolynomialFraction<variable_t> polyFraction = 
      cfd::detail::convert_expression< PolynomialFraction<variable_t> >(normalised);

    return optimised_t(polyFraction);
  }
*/
  ExcafeExpression normalised() const
  {
    return *this;
  }

  std::set<variable_t> getVariables() const
  {
    ExcafeSymbolCollector collector;
    //expr.traverse(collector);

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

/*
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
*/
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
