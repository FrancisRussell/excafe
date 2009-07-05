#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/numeric/math_utilities.hpp>
#include <iostream>

int main(int argc, char** argv)
{
  using cfd::Polynomial; 

  for(std::size_t n=0; n<5; ++n)
  {
    const Polynomial p = cfd::MathUtilities::jacobi(0, 0, n);
    std::cout << "P" << n << "(x): " << p << std::endl;
  }
}
