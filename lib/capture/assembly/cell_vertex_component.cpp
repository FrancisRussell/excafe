#include <ostream>
#include <excafe/capture/assembly/cell_vertex_component.hpp>

namespace excafe
{

namespace detail
{

std::ostream& operator<<(std::ostream& o, const CellVertexComponent& s)
{
  s.write(o);
  return o;
}

}

}
