#ifndef EXCAFE_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_HPP
#define EXCAFE_CAPTURE_ASSEMBLY_SCALAR_PLACEHOLDER_HPP

#include <cstddef>
#include <map>
#include <ostream>
#include <boost/variant.hpp>
#include <boost/operators.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_void.hpp>
#include <boost/blank.hpp>
#include "position_component.hpp"
#include "cell_vertex_component.hpp"
#include "scalar_access.hpp"
#include "basis_coefficient.hpp"
#include "generic_symbol.hpp"
#include "scalar_placeholder_operators.hpp"
#include <excafe/numeric/ginac_expression.hpp>
#include <excafe/numeric/excafe_expression.hpp>

namespace excafe
{

namespace detail
{

namespace
{

struct ScalarPlaceholderExpression
{
  typedef ExcafeExpression<ScalarPlaceholder> type;
};

}

class ScalarPlaceholder : boost::totally_ordered<ScalarPlaceholder>,
  detail::ImmutableArithmetic<ScalarPlaceholderExpression::type, ScalarPlaceholder, ScalarPlaceholder>,
  detail::ImmutableArithmetic<ScalarPlaceholderExpression::type, ScalarPlaceholder, double>,
  detail::ImmutableArithmetic<ScalarPlaceholderExpression::type, double, ScalarPlaceholder>
{
private:
  typedef boost::variant<boost::blank,
                         PositionComponent, 
                         CellVertexComponent, 
                         ScalarAccess, 
                         BasisCoefficient,
                         GenericSymbol> variant_t;
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

  explicit ScalarPlaceholder(const GenericSymbol& c) : value(c)
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
  typename boost::disable_if_c<boost::is_void<typename Visitor::result_type>::value, typename Visitor::result_type>::type 
  apply(Visitor& v) const
  {
    return boost::apply_visitor(v, value);
  }

  template<typename Visitor>
  typename boost::enable_if_c<boost::is_void<typename Visitor::result_type>::value, typename Visitor::result_type>::type 
  apply(Visitor& v) const
  {
    boost::apply_visitor(v, value);
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
