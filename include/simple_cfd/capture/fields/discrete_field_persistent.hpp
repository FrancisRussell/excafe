#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_PERSISTENT_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_PERSISTENT_HPP

#include <string>
#include "discrete_field_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_expr.hpp"
#include "temporal_index_set.hpp"

namespace cfd
{

namespace detail
{

class DiscreteFieldPersistent : public DiscreteFieldExpr
{
private:
  const std::string name;
  const FunctionSpaceExpr::expr_ptr functionSpace;

public:
  DiscreteFieldPersistent(const std::string& _name, const FunctionSpaceExpr::expr_ptr& _functionSpace) :
    name(_name), functionSpace(_functionSpace)
  {
  }

  void accept(DiscreteExprVisitor& visitor)
  {
    visitor.visit(*this);
  }

  virtual TemporalIndexSet getTemporalIndices() const
  {
    return TemporalIndexSet();
  }

  FunctionSpaceExpr::expr_ptr getFunctionSpace() const
  {
    return functionSpace;
  }
};

}

}

#endif
