#include <iostream>
#include <string>
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/numeric/ginac_expression.hpp>
#include <simple_cfd/numeric/excafe_expression.hpp>

template<typename polynomial_t>
void testPolynomial()
{
  std::cout << "5.1: " << polynomial_t(5.1) << std::endl;
  std::cout << "x: " << polynomial_t("x") << std::endl;
  std::cout << "2y: " << polynomial_t(2, "y") << std::endl;
  std::cout << "x^3: " << polynomial_t("x", 3) << std::endl;
  std::cout << "7k^5: " << polynomial_t(7, "k", 5) << std::endl;

  std::cout << std::endl;

  std::cout << "3x^2 - 15: " << polynomial_t(3, "x", 2) - 15 << std::endl;
  std::cout << "2b^2 + 37: " << polynomial_t(2, "b", 2) + 37 << std::endl;
  std::cout << "(a^3 + 36)/7: " << (polynomial_t("a", 3) + 36)/7 << std::endl;
  std::cout << "(b^4 - 17)*11: " << (polynomial_t("b", 4) - 17)*11 << std::endl;

  std::cout << std::endl;

  std::cout << "-(a^2 + 4): " << -(polynomial_t("a", 2) + 4)<< std::endl;

  std::cout << std::endl;

  std::cout << "x^2 + y^2: " << polynomial_t("x", 2) + polynomial_t("y", 2) << std::endl;
  std::cout << "x*y: " << polynomial_t("x")*polynomial_t("y") << std::endl;
  std::cout << "(x^2 + 10) - (y^2 + 5): " << (polynomial_t("x", 2)+10) - (polynomial_t("y", 2) + 5) << std::endl;
  std::cout << "(x-1)*(y-1): " << (polynomial_t("x") - 1)*(polynomial_t("y") - 1) << std::endl;
  std::cout << "(x+y)(x-y): " << (polynomial_t("x") + polynomial_t("y"))*(polynomial_t("x") - polynomial_t("y")) << std::endl;

  std::cout << std::endl;

  const polynomial_t dTest = polynomial_t(0.5, "x", 3) + polynomial_t(4.0, "y", 2);
  std::cout << "(d/dx) 0.5x^3 + 4y^2: " << dTest.derivative("x")  << std::endl;
  std::cout << "(d/dy) 0.5x^3 + 4y^2: " << dTest.derivative("y")  << std::endl;
}

int main(int argc, char** argv)
{
  std::cout << "Testing cfd::Polynomial<std::string>:" << std::endl;
  testPolynomial< cfd::Polynomial<std::string> >();

  std::cout << std::endl << std::endl;

  std::cout << "Testing cfd::GinacExpression<std::string>:" << std::endl;
  testPolynomial< cfd::GinacExpression<std::string> >();

  std::cout << std::endl << std::endl;

  std::cout << "Testing cfd::GinacExpression<std::string>:" << std::endl;
  testPolynomial< cfd::ExcafeExpression<std::string> >();
}
