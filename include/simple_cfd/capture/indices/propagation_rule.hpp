#ifndef SIMPLE_CFD_CAPTURE_INDICES_PROPAGATION_RULE_HPP
#define SIMPLE_CFD_CAPTURE_INDICES_PROPAGATION_RULE_HPP

namespace cfd
{

namespace detail
{

class PropagationRule
{
private:
  DiscreteExpr* const from;
  DiscreteExpr* const to;

public:
  PropagationRule(DiscreteExpr& _from, DiscreteExpr& _to) : from(&_from), to(&_to)
  {
  }

  DiscreteExpr& getFrom() const
  {
    return *from;
  }

  DiscreteExpr& getTo() const
  {
    return *to;
  }

  //FIXME: This ugly pointer-based comparison operator exists so we can use boost::ptr_set.
  bool operator<(const PropagationRule& b) const
  {
    return this < &b;
  }

  virtual ~PropagationRule() {}
};

}

}

#endif
