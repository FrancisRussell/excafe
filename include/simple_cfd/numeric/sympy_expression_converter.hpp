#ifndef SIMPLE_CFD_NUMERIC_SYMPY_EXPRESSION_CONVERTER_HPP
#define SIMPLE_CFD_NUMERIC_SYMPY_EXPRESSION_CONVERTER_HPP

#include <cstddef>
#include <map>
#include <utility>
#include <stack>
#include <sstream>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/object.hpp>
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

template<typename R>
class SymPyExpressionConverter : public NumericExpressionVisitor<typename R::variable_t>
{
private:
  typedef boost::python::object result_type;
  typedef typename R::variable_t variable_t;
  typedef typename NumericExpressionVisitor<variable_t>::value_t value_t;
  std::map<variable_t, std::size_t> variableIDs;
  std::stack<result_type> stack;

  std::size_t getID(const variable_t& var)
  {
    const typename std::map<variable_t, std::size_t>::iterator varIter = variableIDs.find(var);

    if (varIter!=variableIDs.end())
    {
      return varIter->second;
    }
    else
    {
      const std::size_t id = variableIDs.size();
      variableIDs.insert(std::make_pair(var, id));
      return id;
    }
  }

public:
  virtual void visitConstant(const value_t& s)
  
  {
    stack.push(boost::python::make_tuple(CONST, cfd::numeric_cast<double>(s)));
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

  boost::python::dict getNameMap() const
  {
    boost::python::dict result;
    for(typename std::map<variable_t, std::size_t>::const_iterator varIter=variableIDs.begin(); 
      varIter!=variableIDs.end(); ++varIter)
    {
      std::ostringstream nameStream;
      nameStream << varIter->first;
      result[varIter->second] = nameStream.str();
    }
    return result;
  }
};

}
using sympy_bridge::SymPyExpressionConverter;

}
}

#endif
