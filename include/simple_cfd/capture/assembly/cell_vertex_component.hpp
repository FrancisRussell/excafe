#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_CELL_VERTEX_COMPONENT_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_CELL_VERTEX_COMPONENT_HPP

#include <cstddef>
#include <utility>
#include "scalar_placeholder_operators.hpp"

namespace cfd
{

namespace detail
{

class CellVertexComponent : public ScalarPlaceholderOperators<CellVertexComponent>
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
};

}

}

#endif
