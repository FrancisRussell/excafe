#include <ostream>
#include <excafe/capture/assembly/generic_symbol.hpp>

namespace excafe
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
