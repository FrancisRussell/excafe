#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_SCOPING_HPP

namespace cfd
{

namespace detail
{

class DiscreteExprScoping
{
private:
  std::map<TemporalIndexValue*, DiscreteExprScoping> loops;
  std::set<DiscreteExpr*> exprs;

public:
  void addExpressionNode(TemporalIndexSet& indices, DiscreteExpr& expr)
  {
    if (indices.size() == 0)
    {
      exprs.insert(&expr);
    }
    else if (indices.size() == 1)
    {
      // FIXME: Implement me!
    }
  }
};

}

}
#endif
