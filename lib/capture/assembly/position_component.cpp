#include <ostream>
#include <excafe/capture/assembly/position_component.hpp>

namespace excafe
{

namespace detail
{

std::ostream& operator<<(std::ostream& o, const PositionComponent& s)
{
  s.write(o);
  return o;
}

}

}
