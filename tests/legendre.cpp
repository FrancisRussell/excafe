#include <simple_cfd/numeric/excafe_expression.hpp>
#include <simple_cfd/numeric/math_utilities.hpp>
#include <iostream>
#include <string>
#include <iterator>

int main(int argc, char** argv)
{
  using cfd::ExcafeExpression; 

  for(std::size_t n=0; n<5; ++n)
  {
    const std::string x("x");
    const ExcafeExpression<std::string> p = cfd::MathUtilities::jacobi(x, 0, 0, n);
    std::cout << "P" << n << "(x): " << p << std::endl;

    const std::set<double> roots = cfd::MathUtilities::jacobi_roots(0, 0, n);
    std::cout << "P" << n << "(x) roots: ";
    std::copy(roots.begin(), roots.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl << std::endl;
  }
}
