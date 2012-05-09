#ifndef EXCAFE_CAPTURE_ASSEMBLY_POSITION_COMPONENT_HPP
#define EXCAFE_CAPTURE_ASSEMBLY_POSITION_COMPONENT_HPP

#include <cstddef>
#include <ostream>
#include <excafe/numeric/cast.hpp>

namespace excafe
{

namespace detail
{

class PositionComponent
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
    o << "pos[" << excafe::numeric_cast<char>('x'+component) << "]";
  }
};

std::ostream& operator<<(std::ostream& o, const PositionComponent& p);

}

}

#endif
