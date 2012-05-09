#include <cstddef>
#include <iostream>
#include <vector>
#include <string>
#include <excafe/numeric/polynomial.hpp>
#include <excafe/cse/cse_optimiser.hpp>

int main(int argc, char** argv)
{
  typedef excafe::Polynomial<std::string> poly_t;
  std::vector<poly_t> polys;

  const std::string x("x");
  polys.push_back(x + poly_t(-3, x, 3) + poly_t(5, x, 5) + poly_t(-7, x, 7));

  for(std::size_t i=0; i<polys.size(); ++i)
    std::cout << i << ": " << polys[i] << std::endl;

  std::cout << std::endl << "Beginning CSE..." << std::endl;
  excafe::cse::CSEOptimiser<std::string> optimiser(polys.begin(), polys.end());
}
