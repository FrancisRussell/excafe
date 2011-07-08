#ifndef SIMPLE_CFD_NUMERIC_FACTORISER_HPP
#define SIMPLE_CFD_NUMERIC_FACTORISER_HPP

#include <vector>
#include <utility>
#include <simple_cfd/mp/integer.hpp>

namespace cfd
{

class Factoriser
{
public:
  typedef std::pair<mp::Integer, mp::Integer> power_t;
  std::vector<power_t> factor(const mp::Integer& n);

private:
  static const int primes[];
  mp::Integer removeFactor(mp::Integer& value, const mp::Integer& factor) const;
};

}

#endif
