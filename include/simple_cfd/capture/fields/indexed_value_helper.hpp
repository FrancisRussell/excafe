#ifndef SIMPLE_CFD_CAPTURE_FIELDS_INDEXED_VALUE_HELPER_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_INDEXED_VALUE_HELPER_HPP

#include "indexable_value.hpp"
#include "temporal_index_expr.hpp"
#include "discrete_indexed_object.hpp"
#include "discrete_traits.hpp"

namespace cfd
{

namespace detail
{

template<typename discrete_object_tag>
class IndexedValueHelper
{
private:
  typedef IndexableValue<discrete_object_tag> parent_t;
  typedef typename DiscreteTraits<discrete_object_tag>::holder_t holder_t;

  typename parent_t::value_ptr parent;
  const TemporalIndexExpr indexExpr;

public:
  IndexedValueHelper(const typename parent_t::value_ptr& _parent, const TemporalIndexExpr& index) : 
    parent(_parent), indexExpr(index)
  {
  }

  void operator=(const holder_t& expr)
  {
    parent->handleAssignment(indexExpr, expr.getExpr());
  }

  void operator=(const IndexedValueHelper& expr)
  {
    *this = static_cast<holder_t>(expr);
  }

  operator holder_t() const
  {
    return holder_t(new typename DiscreteTraits<discrete_object_tag>::indexed_expr_t(parent, indexExpr));
  }
};

}

}

#endif
