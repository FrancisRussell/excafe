#include <iostream>
#include <string>
#include <boost/python/detail/wrap_python.hpp>
#include <excafe/python_manager.hpp>
#include <excafe/numeric/sympy_expression_converter.hpp>
#include <excafe/numeric/excafe_expression.hpp>
#include <excafe/numeric/sympy_expression.hpp>
#include <excafe/numeric/convert_expression.hpp>

int main(int, char**)
{
  using namespace excafe;

  typedef ExcafeExpression<std::string> poly_t;
  poly_t poly = poly_t(7, "k", 5) + poly_t(4, "x", 2) + abs(poly_t(2, "y"));

  PythonManager::instance().init();

  try
  {
    SymPyExpression<std::string> sympyExpr(poly);
    std::cout << "Excafe expression: " << poly << std::endl;
    std::cout << "After Excafe -> SymPy: " << sympyExpr << std::endl;
    std::cout << "After SymPy -> Excafe: ";
    std::cout << excafe::detail::convert_expression< ExcafeExpression<std::string> >(sympyExpr) << std::endl;

    std::cout << "Integrated w.r.t. x: " << sympyExpr.integrate("x") << std::endl;
    std::cout << "Integrated w.r.t. x (0<x<1): " << sympyExpr.integrate("x", 0, 1) << std::endl;
  }
  catch (boost::python::error_already_set&)
  {
    PyErr_Print();
  }
}
