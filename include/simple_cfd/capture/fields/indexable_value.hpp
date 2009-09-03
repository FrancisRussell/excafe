#ifndef SIMPLE_CFD_CAPTURE_FIELDS_INDEXABLE_VALUE_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_INDEXABLE_VALUE_HPP

#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace detail
{

template<typename discrete_object_tag>
class IndexableValue
{
public:
  typedef boost::shared_ptr<IndexableValue> value_ptr;
  typedef typename DiscreteTraits<discrete_object_tag>::expr_ptr expr_ptr;

  void handleAssignment(const TemporalIndexExpr& indexExpr, const expr_ptr& value)
  {
  }
};

}

}

#endif
