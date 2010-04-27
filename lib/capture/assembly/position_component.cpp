#include <ostream>
#include <simple_cfd/capture/assembly/position_component.hpp>

namespace cfd
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
