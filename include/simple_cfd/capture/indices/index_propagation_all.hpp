#ifndef SIMPLE_CFD_CAPTURE_INDICES_INDEX_PROPAGATION_ALL_HPP
#define SIMPLE_CFD_CAPTURE_INDICES_INDEX_PROPAGATION_ALL_HPP

namespace cfd
{

namespace detail
{

class IndexPropagationAll : public PropagationRule
{
public:
  IndexPropagationAll(DiscreteExpr& from, DiscreteExpr& to) : PropagationRule(from, to)
  {
  }
};

}

}

#endif
