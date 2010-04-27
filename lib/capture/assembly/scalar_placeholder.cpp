#include <ostream>
#include <simple_cfd/capture/assembly/scalar_placeholder.hpp>

namespace cfd
{

namespace detail
{

std::ostream& operator<<(std::ostream& o, const ScalarPlaceholder& s)
{
  s.write(o);
  return o;
}

}

}
