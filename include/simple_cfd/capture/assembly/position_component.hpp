#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_POSITION_COMPONENT_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_POSITION_COMPONENT_HPP

#include <cstddef>
#include <ostream>
#include "scalar_placeholder_operators.hpp"

namespace cfd
{

namespace detail
{

class PositionComponent : public ScalarPlaceholderOperators<PositionComponent>
{
private:
  std::size_t component;

public:
  PositionComponent(const std::size_t _component) : component(_component)
  {
  }

  bool operator==(const PositionComponent& p) const
  {
    return component == p.component;
  }

  bool operator<(const PositionComponent& p) const
  {
    return component < p.component;
  }

  std::size_t getComponent() const
  {
    return component;
  }

  void write(std::ostream& o) const
  {
    o << "pos[" << static_cast<char>('x'+component) << "]";
  }
};

std::ostream& operator<<(std::ostream& o, const PositionComponent& p);

}

}

#endif
