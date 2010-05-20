#include <boost/shared_ptr.hpp>
#include <simple_cfd/capture/fields/temporal_index_value.hpp>
#include <simple_cfd/capture/evaluation/discrete_expr_scoping.hpp>

namespace cfd
{

namespace detail
{

boost::shared_ptr<TemporalIndexValue> DiscreteExprScoping::globalScope(new TemporalIndexValue());

}

}
