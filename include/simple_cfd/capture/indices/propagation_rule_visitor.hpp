#ifndef SIMPLE_CFD_CAPTURE_INDICES_PROPAGATION_RULE_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_INDICES_PROPAGATION_RULE_VISITOR_HPP

#include "indices_fwd.hpp"

namespace cfd
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
