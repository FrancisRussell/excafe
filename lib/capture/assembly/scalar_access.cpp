#include <ostream>
#include <excafe/capture/assembly/scalar_access.hpp>

namespace excafe
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
