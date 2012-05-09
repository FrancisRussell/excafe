#include <ostream>
#include <excafe/capture/assembly/basis_coefficient.hpp>

namespace excafe
{

namespace detail
{

std::ostream& operator<<(std::ostream& o, const BasisCoefficient& c)
{
  c.write(o);
  return o;
}

}

}
