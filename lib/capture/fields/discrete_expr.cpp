#include <excafe/capture/fields/discrete_expr.hpp>
#include <excafe/capture/fields/temporal_index_set.hpp>

namespace excafe
{

namespace detail
{

TemporalIndexSet DiscreteExpr::getTemporalIndices() const
{
  return TemporalIndexSet();
}

TemporalIndexSet DiscreteExpr::getLoopDependencies() const
{
  return TemporalIndexSet();
}

}

}
