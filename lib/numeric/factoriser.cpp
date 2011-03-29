#include <vector>
#include <cln/integer.h>
#include <cln/rational.h>
#include <simple_cfd/numeric/factoriser.hpp>

namespace cfd
{

std::vector<Factoriser::power_t> Factoriser::factor(const cln::cl_I& n)
{
  using namespace cln;

  std::vector< std::pair<cl_I, cl_I> > result;
  cl_I value = n;

  // We explicitly handle the zero case since an empty list of factors is assumed to mean 1.
  if (value == 0)
  {
    result.push_back(power_t(0, 1));
    return result;
  }

  // We generate -1 as a factor for any negative numbers.
  if (value < 0)
  {
    value = -value;
    result.push_back(power_t(-1, 1));
  }

  // Factorise using known primes from table, then revert to using odd numbers.
  std::size_t primeIndex = 0;
  cl_I factor = 0;
  do
  {
    if (primes[primeIndex] != -1)
      factor = primes[primeIndex++];
    else
      factor += 2;

    const cl_I exponent = removeFactor(value, factor);
    if (exponent > 0)
      result.push_back(power_t(factor, exponent));
  }
  while(value > 1);

  return result;
}

cln::cl_I Factoriser::removeFactor(cln::cl_I& value, const cln::cl_I& factor) const
{
  cln::cl_I exponent = 0;
  while (rem(value, factor) == 0)
  {
    value /= factor;
    ++exponent;
  }

  return exponent;
}

}
