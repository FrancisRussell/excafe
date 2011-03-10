#include <simple_cfd/symbolic/basic.hpp>
#include <ostream>

namespace cfd
{

namespace symbolic
{

std::ostream& operator<<(std::ostream& o, const Basic& b)
{
  b.write(o);
  return o;
}

}

}
