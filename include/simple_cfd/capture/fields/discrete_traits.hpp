#ifndef SIMPLE_CFD_CAPTURE_DISCRETE_TRAITS_HPP
#define SIMPLE_CFD_CAPTURE_DISCRETE_TRAITS_HPP

#include "fields_fwd.hpp"
#include "scalar_expr.hpp"
#include "discrete_field_expr.hpp"
#include "operator_expr.hpp"

namespace cfd
{

namespace detail
{

class discrete_scalar_tag {};
class discrete_field_tag {};
class discrete_operator_tag {};

template<typename T>
struct DiscreteTraits
{
};

template<>
struct DiscreteTraits<discrete_scalar_tag>
{
  typedef Scalar holder_t;
  typedef ScalarExpr expr_t;
  typedef ScalarExpr::expr_ptr expr_ptr;
};

template<>
struct DiscreteTraits<discrete_field_tag>
{
  typedef Field holder_t;
  typedef DiscreteFieldExpr expr_t;
  typedef DiscreteFieldExpr::expr_ptr expr_ptr;
};

template<>
struct DiscreteTraits<discrete_operator_tag>
{
  typedef Operator holder_t;
  typedef OperatorExpr expr_t;
  typedef OperatorExpr::expr_ptr expr_ptr;
};

}

}

#endif
