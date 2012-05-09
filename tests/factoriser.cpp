#include <vector>
#include <excafe/numeric/factoriser.hpp>
#include <cln/cln.h>
#include <boost/utility.hpp>
#include <iostream>

void factor(const cln::cl_I& n)
{
  using namespace excafe;

  Factoriser factoriser;
  const std::vector<Factoriser::power_t> result = factoriser.factor(n);

  std::cout << "Factors of " << n << ": ";

  for (std::vector<Factoriser::power_t>::const_iterator 
       iter=result.begin(); iter!=result.end(); 
       ++iter)
  {
    std::cout << iter->first << "^" << iter->second;
    if (boost::next(iter) != result.end())
      std::cout << ", ";
  }
  std::cout << std::endl;
}


int main(int argc, char** argv)
{
  // Factoring 0 should return an explicit 0.
  factor(0);

  // Factoring 1 should return nothing.
  factor(1);

  // Factoring 48 should return 7^2.
  factor(49);

  // Factoring 53 shoudl return itself.
  factor(53);

  // Factoring a negative value.
  factor(153);
  factor(-153);
}
