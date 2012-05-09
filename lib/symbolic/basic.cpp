#include <excafe/symbolic/basic.hpp>
#include <ostream>

namespace excafe
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
