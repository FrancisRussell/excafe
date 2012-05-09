#ifndef EXCAFE_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_VISITOR_HPP
#define EXCAFE_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_VISITOR_HPP

#include "evaluation_fwd.hpp"
#include <excafe/capture/fields/discrete_expr_visitor.hpp>

namespace excafe
{

namespace detail
{

class DiscreteExprScopingVisitor : public DiscreteExprVisitor
{
public:
  virtual void visitBlock(DiscreteExprScoping& scope) = 0;
  virtual void visitLoop(DiscreteExprScoping& scope, TemporalIndexValue* loopIndex) = 0;
};

}

}

#endif
