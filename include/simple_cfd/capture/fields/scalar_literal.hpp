#ifndef SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_LITERAL_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_LITERAL_HPP

#include "scalar_expr.hpp"
#include "discrete_expr_visitor.hpp"

namespace cfd
{

namespace detail
{

class ScalarLiteral : public ScalarExpr
{
private:
  const double value;

public:
  ScalarLiteral(const double _value) : value(_value)
  {
  }

  double getValue() const
  {
    return value;
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
