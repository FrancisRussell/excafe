#include <ostream>
#include <simple_cfd/cse/sop.hpp>

namespace cfd
{

namespace cse
{

std::ostream& operator<<(std::ostream& o, const SOP& sop)
{
  sop.write(o, detail::LiteralWriter());
  return o;
}

}

}
