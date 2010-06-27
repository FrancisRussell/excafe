#include <ostream>
#include <simple_cfd/capture/assembly/scalar_constant.hpp>

namespace cfd
{

namespace detail
{

std::ostream& operator<<(std::ostream& o, const ScalarConstant& c)
{
  c.write(o);
  return o;
}

}

}
