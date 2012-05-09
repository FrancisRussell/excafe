#ifndef EXCAFE_CAPTURE_FIELDS_FIELD_HPP
#define EXCAFE_CAPTURE_FIELDS_FIELD_HPP

#include "fields_fwd.hpp"
#include "discrete_field_expr.hpp"
#include <boost/operators.hpp>

namespace excafe
{

class Field : boost::additive<Field>
{
public:
  typedef detail::DiscreteFieldExpr::expr_ptr expr_ptr;

private:
  expr_ptr expr;

public:
  Field();
  Field(const FunctionSpace& functionSpace);
  Field(detail::DiscreteFieldExpr* const _expr);
  Field(detail::DiscreteFieldExpr::expr_ptr const _expr);
  Field& operator+=(const Field& f);
  Field& operator-=(const Field& f);
  Scalar two_norm() const;
  expr_ptr getExpr() const;
};

}

#endif
