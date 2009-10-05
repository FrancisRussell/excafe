#include <simple_cfd/capture/fields/discrete_expr.hpp>
#include <simple_cfd/capture/fields/temporal_index_set.hpp>

namespace cfd
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
