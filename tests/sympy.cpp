#include <iostream>
#include <string>
#include <simple_cfd/numeric/sympy.hpp>
#include <simple_cfd/numeric/ginac_expression.hpp>

static const char* code =
"def defineSymbols(names):\n"
"  symList = map(lambda (id, name): (id, sympy.Symbol(name)), names.items())\n"
"  return dict(symList)\n"

"def toSymPy(symbolMap, e):\n"
"  if e[0] == SimpleCFD.OperatorType.SYM:\n"
"    return symbolMap[e[1]]\n"
"  elif e[0] == SimpleCFD.OperatorType.CONST:\n"
"    return sympy.Real(e[1])\n"
"  elif e[0] == SimpleCFD.OperatorType.EXP:\n"
"    return toSymPy(symbolMap, e[1][0])**e[1][1]\n"
"  elif e[0] == SimpleCFD.OperatorType.MUL:\n"
"    return reduce((lambda x,y: x*y), map((lambda x: toSymPy(symbolMap, x)), e[1]))\n"
"  elif e[0] == SimpleCFD.OperatorType.ADD:\n"
"    return reduce((lambda x,y: x+y), map((lambda x: toSymPy(symbolMap, x)), e[1]))\n"
"  else:\n"
"    raise Exception('Unknown enum value')\n"

"symbolMap = defineSymbols(symbolNames)\n"
"sympyExpr = toSymPy(symbolMap, expression)\n"
"print sympyExpr\n";

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
