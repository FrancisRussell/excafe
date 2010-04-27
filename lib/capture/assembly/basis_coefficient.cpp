#include <ostream>
#include <simple_cfd/capture/assembly/basis_coefficient.hpp>

namespace cfd
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
