#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_HPP

#include <cstddef>
#include <map>
#include <ostream>
#include <boost/variant.hpp>
#include <boost/operators.hpp>
#include <boost/blank.hpp>
#include "position_component.hpp"
#include "cell_vertex_component.hpp"
#include "scalar_access.hpp"
#include "basis_coefficient.hpp"
#include "scalar_placeholder_operators.hpp"
#include <simple_cfd/numeric/ginac_expression.hpp>

namespace cfd
{

namespace detail
{

namespace
{

struct ScalarPlaceholderExpression
{
  typedef GinacExpression<ScalarPlaceholder> type;
};

}

class ScalarPlaceholder : boost::totally_ordered<ScalarPlaceholder>,
  detail::ImmutableArithmetic<ScalarPlaceholderExpression::type, ScalarPlaceholder, ScalarPlaceholder>,
  detail::ImmutableArithmetic<ScalarPlaceholderExpression::type, ScalarPlaceholder, double>,
  detail::ImmutableArithmetic<ScalarPlaceholderExpression::type, double, ScalarPlaceholder>
{
private:
  typedef boost::variant<boost::blank, PositionComponent, CellVertexComponent, ScalarAccess, BasisCoefficient> variant_t;
  variant_t value;

public:
  typedef ScalarPlaceholderExpression::type expression_t;
  typedef expression_t::optimised_t optimised_expression_t;

  explicit ScalarPlaceholder(const PositionComponent& c) : value(c)
  {
  }

  explicit ScalarPlaceholder(const CellVertexComponent& c) : value(c)
  {
  }

  explicit ScalarPlaceholder(const ScalarAccess& s) : value(s)
  {
  }

  explicit ScalarPlaceholder(const BasisCoefficient& c) : value(c)
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

  void write(std::ostream& o) const
  {
    o << value;
  }
};

std::ostream& operator<<(std::ostream& o, const ScalarPlaceholder& s);

}

}


#endif
