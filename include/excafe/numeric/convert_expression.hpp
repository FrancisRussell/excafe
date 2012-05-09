#ifndef EXCAFE_NUMERIC_EXPRESSION_CONVERTER_HPP
#define EXCAFE_NUMERIC_EXPRESSION_CONVERTER_HPP

#include <stack>
#include <cstddef>
#include <boost/utility/enable_if.hpp>
#include <boost/static_assert.hpp>
#include "cast.hpp"
#include "expression.hpp"
#include "expression_visitor.hpp"
#include "traits.hpp"
#include "symbol_mapper.hpp"
#include "identity_symbol_mapper.hpp"
#include <excafe/exception.hpp>

namespace excafe
{

namespace
{

template<typename R, typename V>
class NumericExpressionConverter : public NumericExpressionVisitor<V>
{
private:
  typedef R result_type;
  typedef typename result_type::variable_t result_variable_t;
  typedef V source_variable_t;

  detail::SymbolMapper<source_variable_t, result_variable_t>& mapper;
  std::stack<result_type> stack;

  template<typename T>
  static typename boost::disable_if_c<T::supports_abs, T>::type performAbs(const T& s)
  {
    CFD_EXCEPTION("Expression conversion result type does not support abs.");
  }

  template<typename T>
  static typename boost::enable_if_c<T::supports_abs, T>::type performAbs(const T& s)
  {
    return abs(s);
  }

public:
  NumericExpressionConverter(detail::SymbolMapper<source_variable_t, result_variable_t>& _mapper) : mapper(_mapper)
  {
  }

  void visitConstant(const typename NumericExpressionVisitor<V>::float_t& s)
  {
    stack.push(result_type(s));
  }

  void visitConstant(const typename NumericExpressionVisitor<V>::integer_t& s)
  {
    stack.push(result_type(s));
  }

  void visitVariable(const source_variable_t& var)
  {
    stack.push(result_type(mapper.getSymbol(var)));
  }

  void visitExponent(const int exponent)
  {
    const result_type val = stack.top();
    stack.pop();
    stack.push(pow(val, exponent));
  }

  void visitAbsoluteValue()
  {
    const result_type val = stack.top();
    stack.pop();
    stack.push(performAbs(val));
  }

  void postSummation(const std::size_t nops)
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

  void postProduct(const std::size_t nops)
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
      CFD_EXCEPTION("Incorrect NumericExpressionVisitor use detected in NumericExpressionConverter.");
    else
      return stack.top();
  }
};

}

namespace detail
{

template<typename result_type, typename source_symbol_type>
result_type convert_expression(const NumericExpression<source_symbol_type>& e, 
                               SymbolMapper<source_symbol_type, typename result_type::variable_t>& mapper)
{
  NumericExpressionConverter<result_type, source_symbol_type> converter(mapper);
  e.accept(converter);
  return converter.getResult();
}

template<typename result_type>
result_type convert_expression(const NumericExpression<typename result_type::variable_t>& e)
{
  detail::IdentitySymbolMapper<typename result_type::variable_t> identityMapper;
  return convert_expression<result_type>(e, identityMapper);
}

}

}

#endif
