#include <simple_cfd/numeric/sympy_bridge.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/exec.hpp>
#include <boost/python/import.hpp>
#include <boost/python/module.hpp>
#include <boost/python/object.hpp>
#include <string>

namespace cfd
{
namespace detail
{
namespace sympy_bridge
{

static const char* initCode =
"def buildSymbolMap(names):\n"
"  symList = map(lambda (id, name): (id, sympy.Symbol(name)), names.items())\n"
"  return dict(symList)\n"

"def symPyZero():\n"
"  return sympy.Real(0.0)\n"

"def toString(o):\n"
"  return str(o)\n"

"def integrate_indefinite(symbolMap, expr, var):\n"
"  return sympy.integrate(expr, symbolMap[var])\n"

"def integrate_definite(symbolMap, expr, var, a, b):\n"
"  return sympy.integrate(expr, (symbolMap[var], a, b))\n"

"def commonToSymPy(symbolMap, e):\n"
"  (op, arg) = e\n"
"  if op == SymPyBridge.OperatorType.SYM:\n"
"    return symbolMap[arg]\n"
"  elif op == SymPyBridge.OperatorType.CONST:\n"
"    return sympy.Real(arg)\n"
"  elif op == SymPyBridge.OperatorType.EXP:\n"
"    return commonToSymPy(symbolMap, arg[0])**arg[1]\n"
"  elif op == SymPyBridge.OperatorType.MUL:\n"
"    return reduce((lambda x,y: x*y), map((lambda x: commonToSymPy(symbolMap, x)), arg))\n"
"  elif op == SymPyBridge.OperatorType.ADD:\n"
"    return reduce((lambda x,y: x+y), map((lambda x: commonToSymPy(symbolMap, x)), arg))\n"
"  else:\n"
"    raise ValueError('Unknown enum value: '+repr(op))\n"

"def symPyToCommon(symbolMap, e):\n"
"  inverseSymbolMap = dict(map((lambda (x,y): (y,x)), symbolMap.items()))\n"
"  return SymPyToCommonConverter().convert(inverseSymbolMap, e)\n"

"class SymPyToCommonConverter(object):\n"
"  def convert(self, symbolMap, expr):\n"
"    exprClass = expr.__class__.__name__\n"
"    convertMethod = '_convert_'+exprClass\n"
"    if hasattr(self, convertMethod):\n"
"      return getattr(self, convertMethod)(symbolMap, expr)\n"
"    else:\n"
"      raise Exception('Cannot convert SymPy type: '+ exprClass)\n"
"\n"
"  def _convert_Add(self, symbolMap, expr):\n"
"    return (SymPyBridge.OperatorType.ADD, map((lambda x: self.convert(symbolMap, x)), expr.args))\n"
"\n"
"  def _convert_Mul(self, symbolMap, expr):\n"
"    return (SymPyBridge.OperatorType.MUL, map((lambda x: self.convert(symbolMap, x)), expr.args))\n"
"\n"
"  def _convert_Real(self, symbolMap, expr):\n"
"    return (SymPyBridge.OperatorType.CONST, float(expr))\n"
"\n"
"  def _convert_Pow(self, symbolMap, expr):\n"
"    if not expr.args[1].is_integer: raise Exception('Cannot convert non-integer exponent')\n"
"    return (SymPyBridge.OperatorType.EXP, (self.convert(symbolMap, expr.args[0]), int(expr.args[1])))\n"
"\n"
"  def _convert_Symbol(self, symbolMap, expr):\n"
"    return (SymPyBridge.OperatorType.SYM, symbolMap[expr])\n";

BOOST_PYTHON_MODULE(SymPyBridge)
{

boost::python::enum_<OperatorType>("OperatorType")
  .value("ADD", ADD)
  .value("MUL", MUL)
  .value("EXP", EXP)
  .value("CONST", CONST)
  .value("SYM", SYM);

}

void init(boost::python::object& global)
{
  initSymPyBridge();
  global["sympy"] = boost::python::import("sympy");
  global["SymPyBridge"] = boost::python::import("SymPyBridge");
  boost::python::exec(initCode, global);
}

}
}
}
