#include <ostream>
#include <simple_cfd/capture/assembly/scalar_access.hpp>

namespace cfd
{

namespace detail
{

std::ostream& operator<<(std::ostream& o, const ScalarAccess& s)
{
  s.write(o);
  return o;
}

}

}
