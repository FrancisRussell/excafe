#ifndef SIMPLE_CFD_NUMERIC_FACTORISER_HPP
#define SIMPLE_CFD_NUMERIC_FACTORISER_HPP

#include <vector>
#include <utility>
#include <cln/integer.h>

namespace cfd
{

class Factoriser
{
public:
  typedef std::pair<cln::cl_I, cln::cl_I> power_t;

  std::vector<power_t> factor(const cln::cl_I& n)
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

    cl_I factor(2);

    while(value != 1)
    {
      cl_I exponent = 0;
      while (rem(value, factor) == 0)
      {
        value /= factor;
        ++exponent;
      }

      if (exponent > 0)
        result.push_back(power_t(factor, exponent));

      ++factor;
    }

    return result;
  }
};

}

#endif
