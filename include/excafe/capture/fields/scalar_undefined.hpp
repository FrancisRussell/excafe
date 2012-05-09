#ifndef EXCAFE_CAPTURE_FIELDS_SCALAR_UNDEFINED_HPP
#define EXCAFE_CAPTURE_FIELDS_SCALAR_UNDEFINED_HPP

#include "discrete_expr_visitor.hpp"
#include "scalar_expr.hpp"

namespace excafe
{

namespace detail
{

class ScalarUndefined : public ScalarExpr
{
public:
  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
