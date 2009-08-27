#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_PERSISTENT_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_PERSISTENT_HPP

#include <string>
#include "discrete_field_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space.hpp"

namespace cfd
{

namespace detail
{

class DiscreteFieldPersistent : public DiscreteFieldExpr
{
private:
  const std::string name;
  const FunctionSpace functionSpace;

public:
  DiscreteFieldPersistent(const std::string& _name, const FunctionSpace& _functionSpace) :
    name(_name), functionSpace(_functionSpace)
  {
  }

  void accept(DiscreteExprVisitor& visitor)
  {
    visitor.visit(*this);
  }
};

}

}

#endif
