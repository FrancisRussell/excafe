#include <simple_cfd/cse/polynomial_index.hpp>

#include <ostream>

namespace cfd
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
