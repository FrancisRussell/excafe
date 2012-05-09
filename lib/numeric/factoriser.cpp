#include <vector>
#include <excafe/numeric/factoriser.hpp>
#include <excafe/mp/integer.hpp>

namespace excafe
{

std::vector<Factoriser::power_t> Factoriser::factor(const mp::Integer& n)
{
  using namespace mp;

  std::vector< std::pair<Integer, Integer> > result;
  Integer value = n;

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

  // Compute upper bound on factors
  Integer sqrtFloor = isqrt(value);

  // Factorise using known primes from table, then revert to using odd numbers.
  std::size_t primeIndex = 0;
  Integer factor = 0;
  do
  {
    if (primes[primeIndex] != -1)
      factor = primes[primeIndex++];
    else
      factor += 2;

    const Integer exponent = removeFactor(value, factor);
    if (exponent > 0)
    {
      result.push_back(power_t(factor, exponent));
      sqrtFloor = isqrt(value);
    }
  }
  while(value != 1 && factor <= sqrtFloor);

  // value may still contain the final prime factor
  if (value != 1)
    result.push_back(power_t(value, 1));

  return result;
}

mp::Integer Factoriser::removeFactor(mp::Integer& value, const mp::Integer& factor) const
{
  mp::Integer exponent = 0;
  while ((value%factor) == 0)
  {
    value /= factor;
    ++exponent;
  }

  return exponent;
}

}
