#include <iostream>
#include <string>
#include <boost/python/detail/wrap_python.hpp>
#include <simple_cfd/python_manager.hpp>
#include <simple_cfd/numeric/sympy_expression_converter.hpp>
#include <simple_cfd/numeric/ginac_expression.hpp>
#include <simple_cfd/numeric/sympy_expression.hpp>
#include <simple_cfd/numeric/convert_expression.hpp>

int main(int, char**)
{
  using namespace cfd;

  typedef GinacExpression<std::string> poly_t;
  poly_t poly = poly_t(7, "k", 5) + poly_t(4, "x", 2);

  PythonManager::instance().init();

  try
  {
    SymPyExpression<std::string> sympyExpr(poly);
    std::cout << "GiNaC expression: " << poly << std::endl;
    std::cout << "After GiNaC -> SymPy: " << sympyExpr << std::endl;
    std::cout << "After SymPy -> GiNaC: ";
    std::cout << cfd::detail::convert_expression< GinacExpression<std::string> >(sympyExpr) << std::endl;

    std::cout << "Integrated w.r.t. x: " << sympyExpr.integrate("x") << std::endl;
    std::cout << "Integrated w.r.t. x (0<x<1): " << sympyExpr.integrate("x", 0, 1) << std::endl;
  }
  catch (boost::python::error_already_set&)
  {
    PyErr_Print();
  }
}
