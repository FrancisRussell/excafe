#include <ostream>
#include <simple_cfd/capture/assembly/cell_vertex_component.hpp>

namespace cfd
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
