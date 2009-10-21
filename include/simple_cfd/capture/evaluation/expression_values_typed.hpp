#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_EXPRESSION_VALUES_TYPED_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_EXPRESSION_VALUES_TYPED_HPP

#include <cstddef>
#include <map>
#include <utility>
#include <simple_cfd/discrete_value_traits.hpp>

namespace cfd
{

namespace detail
{

template<typename discrete_object_tag, std::size_t D>
class ExpressionValuesTyped
{
private:
  static const std::size_t dimension = D;
  typedef typename DiscreteTraits<discrete_object_tag>::expr_t expr_t;
  typedef typename DiscreteTraits<discrete_object_tag>::indexed_expr_t indexed_expr_t;
  typedef typename DiscreteValueTraits<discrete_object_tag, dimension>::value_t value_t;

  std::map<expr_t*, value_t> nonIndexed;
  std::map<std::pair<indexed_expr_t*, int>, value_t> indexed;
};

}

}

#endif
