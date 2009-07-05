#include <iostream>
#include <simple_cfd/numeric/polynomial.hpp>


int main(int argc, char** argv)
{
  using cfd::Polynomial;

  std::cout << "5.1: " << Polynomial(5.1) << std::endl;
  std::cout << "x: " << Polynomial("x") << std::endl;
  std::cout << "2y: " << Polynomial(2, "y") << std::endl;
  std::cout << "x^3: " << Polynomial("x", 3) << std::endl;
  std::cout << "7k^5: " << Polynomial(7, "k", 5) << std::endl;

  std::cout << std::endl;

  std::cout << "3x^2 - 15: " << Polynomial(3, "x", 2) - 15 << std::endl;
  std::cout << "2b^2 + 37: " << Polynomial(2, "b", 2) + 37 << std::endl;
  std::cout << "(a^3 + 36)/7: " << (Polynomial("a", 3) + 36)/7 << std::endl;
  std::cout << "(b^4 - 17)*11: " << (Polynomial("b", 4) - 17)*11 << std::endl;

  std::cout << std::endl;

  std::cout << "-(a^2 + 4): " << -(Polynomial("a", 2) + 4)<< std::endl;

  std::cout << std::endl;

  std::cout << "x^2 + y^2: " << Polynomial("x", 2) + Polynomial("y", 2) << std::endl;
  std::cout << "x*y: " << Polynomial("x")*Polynomial("y") << std::endl;
  std::cout << "(x^2 + 10) - (y^2 + 5): " << (Polynomial("x", 2)+10) - (Polynomial("y", 2) + 5) << std::endl;
  std::cout << "(x-1)*(y-1): " << (Polynomial("x") - 1)*(Polynomial("y") - 1) << std::endl;
  std::cout << "(x+y)(x-y): " << (Polynomial("x") + Polynomial("y"))*(Polynomial("x") - Polynomial("y")) << std::endl;

  std::cout << std::endl;


  const Polynomial dTest = Polynomial(0.5, "x", 3) + Polynomial(4.0, "y", 2);
  std::cout << "(d/dx) 0.5x^3 + 4y^2: " << dTest.derivative("x")  << std::endl;
  std::cout << "(d/dy) 0.5x^3 + 4y^2: " << dTest.derivative("y")  << std::endl;

}
