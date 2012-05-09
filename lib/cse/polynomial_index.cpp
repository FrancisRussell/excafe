#include <excafe/cse/polynomial_index.hpp>

#include <ostream>

namespace excafe
{

namespace cse
{

std::ostream& operator<<(std::ostream& o, const PolynomialIndex& i)
{
  i.write(o);
  return o;
}

}

}
