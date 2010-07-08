#include <iostream>
#include <string>
#include <boost/python/detail/wrap_python.hpp>
#include <simple_cfd/python_manager.hpp>
#include <simple_cfd/numeric/sympy_expression_converter.hpp>
#include <simple_cfd/numeric/ginac_expression.hpp>

static const char* code =
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

  PythonManager::instance().init();

  try
  {
    cfd::detail::SymPyExpressionConverter< GinacExpression<std::string> > converter;
    poly.accept(converter);

    dict local;
    local["symbolNames"] = converter.getNameMap();
    local["expression"] = converter.getResult();
    object result = PythonManager::instance().execute(code, local);
  }
  catch (error_already_set)
  {
    PyErr_Print();
  }
}
