#include <iostream>
#include <string>
#include <simple_cfd/numeric/sympy.hpp>
#include <simple_cfd/numeric/ginac_expression.hpp>

static const char* code =
"def defineSymbols(names):\n"
"  symList = map(lambda (id, name): (id, sympy.Symbol(name)), names.items())\n"
"  return dict(symList)\n"

"def commonToSymPy(symbolMap, e):\n"
"  op = e[0]\n"
"  arg = e[1]\n"
"  if op == SimpleCFD.OperatorType.SYM:\n"
"    return symbolMap[arg]\n"
"  elif op == SimpleCFD.OperatorType.CONST:\n"
"    return sympy.Real(arg)\n"
"  elif op == SimpleCFD.OperatorType.EXP:\n"
"    return commonToSymPy(symbolMap, arg[0])**arg[1]\n"
"  elif op == SimpleCFD.OperatorType.MUL:\n"
"    return reduce((lambda x,y: x*y), map((lambda x: commonToSymPy(symbolMap, x)), arg))\n"
"  elif op == SimpleCFD.OperatorType.ADD:\n"
"    return reduce((lambda x,y: x+y), map((lambda x: commonToSymPy(symbolMap, x)), arg))\n"
"  else:\n"
"    raise Exception('Unknown enum value')\n"

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
"    return (SimpleCFD.OperatorType.ADD, map((lambda x: self.convert(symbolMap, x)), expr.args))\n"
"\n"
"  def _convert_Mul(self, symbolMap, expr):\n"
"    return (SimpleCFD.OperatorType.MUL, map((lambda x: self.convert(symbolMap, x)), expr.args))\n"
"\n"
"  def _convert_Real(self, symbolMap, expr):\n"
"    return (SimpleCFD.OperatorType.CONST, float(expr))\n"
"\n"
"  def _convert_Pow(self, symbolMap, expr):\n"
"    if not expr.args[1].is_integer: raise Exception('Cannot convert non-integer exponent')\n"
"    return (SimpleCFD.OperatorType.EXP, (self.convert(symbolMap, expr.args[0]), expr.args[1]))\n"
"\n"
"  def _convert_Symbol(self, symbolMap, expr):\n"
"    return (SimpleCFD.OperatorType.SYM, symbolMap[expr])\n"

"symbolMap = defineSymbols(symbolNames)\n"
"sympyExpr = commonToSymPy(symbolMap, expression)\n"
"commonExpr = symPyToCommon(symbolMap, sympyExpr)\n"
"print sympyExpr\n"
"print commonExpr\n";

int main(int, char**)
{
  using namespace cfd;
  using namespace boost::python;

  typedef GinacExpression<std::string> poly_t;
  poly_t poly = poly_t(7, "k", 5) + poly_t(4, "x", 2);
  Py_Initialize();

  try
  {
    cfd::sympy::initSimpleCFD();

    sympy::SymPyExpressionConverter< GinacExpression<std::string> > converter;
    poly.accept(converter);

    object main = import("__main__");
    object global(main.attr("__dict__"));
    global["sympy"] = import("sympy");
    global["SimpleCFD"] = import("SimpleCFD");
    global["symbolNames"] = converter.getNameMap();
    global["expression"] = converter.getResult();
    object result = exec(code, global, global);
  }
  catch (error_already_set)
  {
    PyErr_Print();
  }
}
