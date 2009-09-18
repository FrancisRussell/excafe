#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_TWO_NORM_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_FIELD_TWO_NORM_HPP

#include "scalar_expr.hpp"
#include "discrete_field_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "temporal_index_set.hpp"

namespace cfd
{

namespace detail
{

class DiscreteFieldTwoNorm : public ScalarExpr
{
private:
  DiscreteFieldExpr::expr_ptr field;

public:
  DiscreteFieldTwoNorm(const DiscreteFieldExpr::expr_ptr& f) : field(f)
  {
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual TemporalIndexSet getTemporalIndices() const
  {
    return field->getTemporalIndices();
  }

  
  DiscreteFieldExpr& getField() const
  {
    return *field;
  }
};

}

}

#endif
