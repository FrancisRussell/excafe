#include <cstddef>
#include <iostream>
#include <vector>
#include <string>
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/cse/cse_optimiser.hpp>

int main(int argc, char** argv)
{
  typedef cfd::Polynomial<std::string> poly_t;
  std::vector<poly_t> polys;

  polys.push_back(poly_t("x") + poly_t(-3, "x", 3) + poly_t(5, "x", 5) + poly_t(-7, "x", 7));

  for(std::size_t i=0; i<polys.size(); ++i)
  {
    std::cout << i << ": " << polys[i] << std::endl;
  }

  std::cout << std::endl << "Beginning CSE..." << std::endl;
  cfd::cse::CSEOptimiser<std::string> optimiser(polys.begin(), polys.end());
}
