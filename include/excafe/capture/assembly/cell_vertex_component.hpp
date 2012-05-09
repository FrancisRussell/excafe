#ifndef EXCAFE_CAPTURE_ASSEMBLY_CELL_VERTEX_COMPONENT_HPP
#define EXCAFE_CAPTURE_ASSEMBLY_CELL_VERTEX_COMPONENT_HPP

#include <cstddef>
#include <ostream>
#include <utility>
#include "assembly_fwd.hpp"
#include <excafe/numeric/cast.hpp>

namespace excafe
{

namespace detail
{

class CellVertexComponent
{
private:
  std::size_t vertex;
  std::size_t component;

public:
  CellVertexComponent(const std::size_t _vertex, const std::size_t _component) : 
    vertex(_vertex), component(_component)
  {
  }

  bool operator==(const CellVertexComponent& p) const
  {
    return vertex == p.vertex &&
           component == p.component;
  }

  bool operator<(const CellVertexComponent& p) const
  {
    return std::make_pair(vertex, component) < 
           std::make_pair(p.vertex, p.component);
  }

  std::size_t getVertexID() const
  {
    return vertex;
  }

  std::size_t getComponent() const
  {
    return component;
  }

  void write(std::ostream& o) const
  {
    o << "cell[" << vertex << "][" << excafe::numeric_cast<char>('x'+component) << "]";
  }
};

std::ostream& operator<<(std::ostream& o, const CellVertexComponent& c);

}

}

#endif
