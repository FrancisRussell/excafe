#ifndef SIMPLE_CFD_NUMERIC_SYMPY_EXPRESSION_CONVERTER_HPP
#define SIMPLE_CFD_NUMERIC_SYMPY_EXPRESSION_CONVERTER_HPP

#include <cstddef>
#include <utility>
#include <stack>
#include <sstream>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/object.hpp>
#include <boost/bimap.hpp>
#include "expression_visitor.hpp"
#include "cast.hpp"
#include "sympy_bridge.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{
namespace detail
{
namespace sympy_bridge
{

template<typename V>
class SymPyExpressionConverter : public NumericExpressionVisitor<V>
{
private:
  typedef V variable_t;
  typedef boost::python::object result_type;
  typedef typename NumericExpressionVisitor<variable_t>::float_t   float_t;
  typedef typename NumericExpressionVisitor<variable_t>::integer_t integer_t;
  typedef boost::bimap<variable_t, std::size_t> variable_map_t;
  variable_map_t variableIDs;
  std::stack<result_type> stack;

  std::size_t getID(const variable_t& var)
  {
    const typename variable_map_t::left_iterator varIter = variableIDs.left.find(var);

    if (varIter!=variableIDs.left.end())
    {
      return varIter->second;
    }
    else
    {
      const std::size_t id = variableIDs.size();
      variableIDs.left.insert(std::make_pair(var, id));
      return id;
    }
  }

public:
  virtual void visitConstant(const float_t& s)
  {
    stack.push(boost::python::make_tuple(FLOAT, cfd::numeric_cast<double>(s)));
  }

  virtual void visitConstant(const integer_t& s)
  {
    stack.push(boost::python::make_tuple(INTEGER, cfd::numeric_cast<long>(s)));
  }

  virtual void visitVariable(const variable_t& var)
  {
    stack.push(boost::python::make_tuple(SYM, getID(var)));
  }

  virtual void visitExponent(const int exponent)
  {
    const result_type val = stack.top();
    stack.pop();

    stack.push(boost::python::make_tuple(EXP, boost::python::make_tuple(val, exponent)));
  }

  virtual void visitAbsoluteValue()
  {
    const result_type val = stack.top();
    stack.pop();

    stack.push(boost::python::make_tuple(ABS, val));
  }

  virtual void postSummation(const std::size_t nops)
  {
    if (stack.size() < nops)
      CFD_EXCEPTION("Too few operands on stack for sum in SymPyExpressionConverter.");

    boost::python::list operands;
    for(std::size_t i=0; i<nops; ++i)
    {
      operands.append(stack.top());
      stack.pop();
    }
    
    stack.push(boost::python::make_tuple(ADD, operands));
  }

  virtual void postProduct(const std::size_t nops)
  {
    if (stack.size() < nops)
      CFD_EXCEPTION("Too few operands on stack for product in SymPyExpressionConverter.");

    boost::python::list operands;
    for(std::size_t i=0; i<nops; ++i)
    {
      operands.append(stack.top());
      stack.pop();
    }

    stack.push(boost::python::make_tuple(MUL, operands));
  }

  result_type getResult() const
  {
    if (stack.size() != 1)
    {
      CFD_EXCEPTION("Incorrect NumericExpressionVisitor use detected in SymPyExpressionConverter.");
    }
    else
    {
      return stack.top();
    }
  }

  boost::bimap<variable_t, std::size_t> getVariableMap() const
  {
    return variableIDs;
  }

  boost::python::dict getNameMap() const
  {
    boost::python::dict result;
    for(typename variable_map_t::const_iterator varIter=variableIDs.begin(); 
      varIter!=variableIDs.end(); ++varIter)
    {
      std::ostringstream nameStream;
      nameStream << varIter->left;
      result[varIter->right] = nameStream.str();
    }
    return result;
  }
};
}
using sympy_bridge::SymPyExpressionConverter;

}
}

#endif
