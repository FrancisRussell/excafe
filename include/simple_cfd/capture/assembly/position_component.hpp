#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_POSITION_COMPONENT_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_POSITION_COMPONENT_HPP

#include <cstddef>
#include <ostream>
#include "scalar_placeholder_operators.hpp"
#include <simple_cfd/numeric/cast.hpp>

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
    o << "pos[" << cfd::numeric_cast<char>('x'+component) << "]";
  }
};

std::ostream& operator<<(std::ostream& o, const PositionComponent& p);

}

}

#endif
