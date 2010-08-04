#ifndef SIMPLE_CFD_NUMERIC_EXPRESSION_CONVERTER_HPP
#define SIMPLE_CFD_NUMERIC_EXPRESSION_CONVERTER_HPP

#include <stack>
#include <cstddef>
#include "cast.hpp"
#include "expression.hpp"
#include "expression_visitor.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace
{

template<typename R>
class NumericExpressionConverter : public NumericExpressionVisitor<typename R::variable_t>
{
private:
  typedef R result_type;
  typedef typename result_type::variable_t variable_t;
  typedef typename NumericExpressionVisitor<variable_t>::value_t value_t;
  std::stack<result_type> stack;

public:
  virtual void visitConstant(const value_t& s)
  {
    stack.push(result_type(cfd::numeric_cast<typename result_type::value_type>(s)));
  }

  virtual void visitVariable(const variable_t& var)
  {
    stack.push(result_type(var));
  }

  virtual void visitExponent(const int exponent)
  {
    const result_type val = stack.top();
    stack.pop();

    result_type result(1.0);

    for(std::size_t i=0; i<std::abs(exponent); ++i)
      result *= val;

    stack.push(exponent>0 ? result : 1.0/result);
  }

  virtual void postSummation(const std::size_t nops)
  {
    if (stack.size() < nops)
      CFD_EXCEPTION("Too few operands on stack for sum in NumericExpressionConverter.");

    result_type result = 0.0;
    for(std::size_t i=0; i<nops; ++i)
    {
      result += stack.top();
      stack.pop();
    }
    
    stack.push(result);
  }

  virtual void postProduct(const std::size_t nops)
  {
    if (stack.size() < nops)
      CFD_EXCEPTION("Too few operands on stack for product in NumericExpressionConverter.");

    result_type result = 1.0;
    for(std::size_t i=0; i<nops; ++i)
    {
      result *= stack.top();
      stack.pop();
    }

    stack.push(result);
  }

  result_type getResult() const
  {
    if (stack.size() != 1)
    {
      CFD_EXCEPTION("Incorrect NumericExpressionVisitor use detected in NumericExpressionConverter.");
    }
    else
    {
      return stack.top();
    }
  }
};

}

namespace detail
{

template<typename result_type>
result_type convert_expression(const NumericExpression<typename result_type::variable_t>& e)
{
  NumericExpressionConverter<result_type> converter;
  e.accept(converter);
  return converter.getResult();
}

}

}

#endif