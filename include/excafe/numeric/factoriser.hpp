#ifndef EXCAFE_NUMERIC_FACTORISER_HPP
#define EXCAFE_NUMERIC_FACTORISER_HPP

#include <vector>
#include <utility>
#include <map>
#include <excafe/mp/integer.hpp>

namespace excafe
{

class Factoriser
{
public:
  typedef std::pair<mp::Integer, mp::Integer> power_t;
  std::vector<power_t> factor(const mp::Integer& n);
  std::map< mp::Integer, std::vector<power_t> > cache;

private:
  static const int primes[];
  mp::Integer removeFactor(mp::Integer& value, const mp::Integer& factor) const;
};

}

#endif
