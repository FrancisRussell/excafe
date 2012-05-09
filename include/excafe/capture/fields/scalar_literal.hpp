#ifndef EXCAFE_CAPTURE_FIELDS_SCALAR_LITERAL_HPP
#define EXCAFE_CAPTURE_FIELDS_SCALAR_LITERAL_HPP

#include "scalar_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "temporal_index_set.hpp"
#include <excafe/capture/indices/propagation_rules.hpp>

namespace excafe
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

  PropagationRules getPropagationRules()
  {
    return PropagationRules();
  }

  virtual std::set<DiscreteExpr*> getDependencies() const 
  {
    return std::set<DiscreteExpr*>();
  }
};

}

}

#endif
