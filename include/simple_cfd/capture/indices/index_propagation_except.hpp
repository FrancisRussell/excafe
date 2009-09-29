#ifndef SIMPLE_CFD_CAPTURE_INDICES_INDEX_PROPAGATION_EXCEPT_HPP
#define SIMPLE_CFD_CAPTURE_INDICES_INDEX_PROPAGATION_EXCEPT_HPP

namespace cfd
{

namespace detail
{

class IndexPropagationExcept : public PropagationRule
{
private:
  TemporalIndexValue* const index;

public:
  IndexPropagationExcept(DiscreteExpr& from, DiscreteExpr& to, TemporalIndexValue& i) : 
    PropagationRule(from, to), index(&i)
  {
  }
};

}

}

#endif
