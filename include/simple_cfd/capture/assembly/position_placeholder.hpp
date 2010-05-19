#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_POSITION_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_POSITION_PLACEHOLDER_HPP

#include <cstddef>
#include "scalar_placeholder.hpp"
#include "position_component.hpp"

namespace cfd
{

namespace detail
{

class PositionPlaceholder
{
public:
  ScalarPlaceholder operator[](const std::size_t d) const
  {
    return ScalarPlaceholder(PositionComponent(d));
  }
};

}

}

#endif
