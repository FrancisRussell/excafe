#ifndef SIMPLE_CFD_NUMERIC_SYMPY_EXPRESSION_HPP
#define SIMPLE_CFD_NUMERIC_SYMPY_EXPRESSION_HPP

#include <boost/bimap.hpp>
#include <boost/python/object.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <simple_cfd/numeric/expression.hpp>
#include <simple_cfd/numeric/expression_visitor.hpp>
#include <simple_cfd/numeric/sympy_bridge.hpp>
#include <simple_cfd/numeric/cast.hpp>
#include <simple_cfd/exception.hpp>

namespace cfd
{

template<typename V>
class SymPyExpression : public NumericExpression<V>
{
public:
  typedef V variable_t;

private:
  boost::bimap<variable_t, std::size_t> variableMap;
  boost::python::object symbolMap;
  boost::python::object sympyExpression;

  std::size_t visitList(NumericExpressionVisitor<V>& v, const boost::python::list& list) const
  {
    const std::size_t length = boost::python::len(list);
    for(std::size_t i=0; i<length; ++i)
      visitTuple(v, boost::python::extract<boost::python::tuple>(list[i]));
    return length;
  }

  void visitTuple(NumericExpressionVisitor<V>& v, const boost::python::tuple& tuple) const
  {
    using namespace detail::sympy_bridge;

    const OperatorType opType = boost::python::extract<OperatorType>(tuple[0]);
    const boost::python::object operand = tuple[1];

    if (opType == ADD)
    {
      const std::size_t nops = visitList(v, boost::python::extract<boost::python::list>(operand));
      v.postSummation(nops);
    }
    else if (opType == MUL)
    {
      const std::size_t nops = visitList(v, boost::python::extract<boost::python::list>(operand));
      v.postProduct(nops);
    }
    else if (opType == EXP)
    {
      const boost::python::tuple expTuple = boost::python::extract<boost::python::tuple>(operand);
      visitTuple(v, boost::python::extract<boost::python::tuple>(expTuple[0]));
      v.visitExponent(boost::python::extract<int>(expTuple[1]));
    }
    else if (opType == CONST)
    {
      const double constant = boost::python::extract<double>(operand);
      v.visitConstant(cfd::numeric_cast<typename NumericExpressionVisitor<V>::value_t>(constant));
    }
    else if (opType == SYM)
    {
      const std::size_t symID = boost::python::extract<std::size_t>(operand);
      typename boost::bimap<variable_t, std::size_t>::right_const_iterator iter = variableMap.right.find(symID);
      if (iter != variableMap.right.end())
      {
        v.visitVariable(iter->second);
      }
      else
      {
        CFD_EXCEPTION("Unknown symbol identifier passed back from SymPy.");
      }
    }
    else
    {
      CFD_EXCEPTION("Unknown operation type when parsing marshalled SymPy expression.");
    }
  }

public:
  SymPyExpression()
  {
    const boost::python::object buildZero = PythonManager::instance().getObject("symPyZero");
    sympyExpression = buildZero();
  }

  SymPyExpression(const NumericExpression<variable_t>& source)
  {
    detail::SymPyExpressionConverter<variable_t> converter;
    source.accept(converter);

    const boost::python::object buildSymbolMap = PythonManager::instance().getObject("buildSymbolMap");
    const boost::python::object commonToSymPy = PythonManager::instance().getObject("commonToSymPy");
    symbolMap = buildSymbolMap(converter.getNameMap());
    sympyExpression = commonToSymPy(symbolMap, converter.getResult());
    variableMap = converter.getVariableMap();
  }

  void write(std::ostream& o) const
  {
    const boost::python::object str = PythonManager::instance().getObject("toString");
    const char* desc = boost::python::extract<const char*>(str(sympyExpression));
    o << desc;
  }
  
  virtual void accept(NumericExpressionVisitor<V>& v) const
  {
    const boost::python::object symPyToCommon = PythonManager::instance().getObject("symPyToCommon");
    const boost::python::object expression = symPyToCommon(symbolMap, sympyExpression);
    visitTuple(v, boost::python::extract<boost::python::tuple>(expression));
  }
};

template<typename V>
std::ostream& operator<<(std::ostream& o, const SymPyExpression<V>& e)
{
  e.write(o);
  return o;
}

}

#endif
