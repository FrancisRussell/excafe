#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_HPP

#include <cstddef>
#include <map>
#include <boost/variant.hpp>
#include <boost/operators.hpp>
#include "position_component.hpp"
#include "cell_vertex_component.hpp"
#include "scalar_access.hpp"
#include "basis_coefficient.hpp"

namespace cfd
{

namespace detail
{

class ScalarPlaceholder : boost::equality_comparable<ScalarPlaceholder>
{
private:
  typedef boost::variant<PositionComponent, CellVertexComponent, ScalarAccess, BasisCoefficient> variant_t;
  variant_t value;

public:
  ScalarPlaceholder(const PositionComponent& c) : value(c)
  {
  }

  ScalarPlaceholder(const CellVertexComponent& c) : value(c)
  {
  }

  ScalarPlaceholder(const ScalarAccess& s) : value(s)
  {
  }

  ScalarPlaceholder(const BasisCoefficient& c) : value(c)
  {
  }

  bool operator==(const ScalarPlaceholder& s) const
  {
    return value == s.value;
  }

  bool operator<(const ScalarPlaceholder& s) const
  {
    return value < s.value;
  }

  template<typename Visitor>
  typename Visitor::result_type apply(const Visitor& v) const
  {
    return boost::apply_visitor(v, value);
  }
};

}

}

#endif
