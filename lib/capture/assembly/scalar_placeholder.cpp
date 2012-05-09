#include <ostream>
#include <excafe/capture/assembly/scalar_placeholder.hpp>

namespace excafe
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
