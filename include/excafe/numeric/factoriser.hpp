#ifndef EXCAFE_NUMERIC_FACTORISER_HPP
#define EXCAFE_NUMERIC_FACTORISER_HPP

#include <vector>
#include <utility>
#include <cln/integer.h>

namespace excafe
{

class Factoriser
{
public:
  typedef std::pair<cln::cl_I, cln::cl_I> power_t;
  std::vector<power_t> factor(const cln::cl_I& n);

private:
  static const int primes[];
  cln::cl_I removeFactor(cln::cl_I& value, const cln::cl_I& factor) const;
};

}

#endif
