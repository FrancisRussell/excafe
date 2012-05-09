#ifndef EXCAFE_CAPTURE_FIELDS_DISCRETE_TRAITS_HPP
#define EXCAFE_CAPTURE_FIELDS_DISCRETE_TRAITS_HPP

#include "fields_fwd.hpp"
#include "scalar_expr.hpp"
#include "discrete_field_expr.hpp"
#include "operator_expr.hpp"
#include "scalar_undefined.hpp"
#include "discrete_field_undefined.hpp"
#include "operator_undefined.hpp"

namespace excafe
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
  typedef ScalarExpr::weak_expr_ptr weak_expr_ptr;
  typedef DiscreteIndexedScalar indexed_expr_t;  
  typedef ScalarUndefined undefined_t;  
  typedef IndexableValue<discrete_scalar_tag> indexable_t;
};

template<>
struct DiscreteTraits<discrete_field_tag>
{
  typedef Field holder_t;
  typedef DiscreteFieldExpr expr_t;
  typedef DiscreteFieldExpr::expr_ptr expr_ptr;
  typedef DiscreteFieldExpr::weak_expr_ptr weak_expr_ptr;
  typedef DiscreteIndexedField indexed_expr_t;  
  typedef DiscreteFieldUndefined undefined_t;  
  typedef IndexableValue<discrete_field_tag> indexable_t;
};

template<>
struct DiscreteTraits<discrete_operator_tag>
{
  typedef Operator holder_t;
  typedef OperatorExpr expr_t;
  typedef OperatorExpr::expr_ptr expr_ptr;
  typedef OperatorExpr::weak_expr_ptr weak_expr_ptr;
  typedef DiscreteIndexedOperator indexed_expr_t;  
  typedef OperatorUndefined undefined_t;  
  typedef IndexableValue<discrete_operator_tag> indexable_t;
};

}

}

#endif
