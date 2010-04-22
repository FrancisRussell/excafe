#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_POSITION_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_POSITION_PLACEHOLDER_HPP

#include <cstddef>
#include "position_component.hpp"

namespace cfd
{

namespace detail
{

class PositionPlaceholder
{
public:
  PositionComponent operator[](const std::size_t d) const
  {
    return PositionComponent(d);
  }
};

}

}

#endif
