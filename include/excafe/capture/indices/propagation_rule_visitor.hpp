#ifndef EXCAFE_CAPTURE_INDICES_PROPAGATION_RULE_VISITOR_HPP
#define EXCAFE_CAPTURE_INDICES_PROPAGATION_RULE_VISITOR_HPP

#include "indices_fwd.hpp"

namespace excafe
{

namespace detail
{

class PropagationRuleVisitor
{
public:
  virtual void visit(IndexPropagationAll& i) = 0;
  virtual void visit(IndexPropagationExcept& i) = 0;
  virtual ~PropagationRuleVisitor() {}
};

}

}

#endif
