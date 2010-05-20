#ifndef SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_UNDEFINED_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_UNDEFINED_HPP

#include "discrete_expr_visitor.hpp"
#include "scalar_expr.hpp"

namespace cfd
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
