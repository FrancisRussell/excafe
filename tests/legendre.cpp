#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/numeric/math_utilities.hpp>
#include <iostream>
#include <string>
#include <iterator>

int main(int argc, char** argv)
{
  using cfd::Polynomial; 

  for(std::size_t n=0; n<5; ++n)
  {
    const Polynomial<std::string> p = cfd::MathUtilities::jacobi(0, 0, n);
    std::cout << "P" << n << "(x): " << p << std::endl;

    const std::set<double> roots = cfd::MathUtilities::jacobi_roots(0, 0, n);
    std::cout << "P" << n << "(x) roots: ";
    std::copy(roots.begin(), roots.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl << std::endl;
  }
}
