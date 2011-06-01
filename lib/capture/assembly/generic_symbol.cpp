#include <ostream>
#include <simple_cfd/capture/assembly/generic_symbol.hpp>

namespace cfd
{

namespace detail
{

std::ostream& operator<<(std::ostream& o, const GenericSymbol& s)
{
  s.write(o);
  return o;
}

}

}
