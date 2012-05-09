#ifndef EXCAFE_CAPTURE_FIELDS_INDEXED_VALUE_HELPER_HPP
#define EXCAFE_CAPTURE_FIELDS_INDEXED_VALUE_HELPER_HPP

#include "indexable_value.hpp"
#include "temporal_index_expr.hpp"
#include "discrete_indexed_object.hpp"
#include "discrete_traits.hpp"

namespace excafe
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
    typename DiscreteTraits<discrete_object_tag>::expr_ptr expr = parent->getIndexedExpr(indexExpr);

    // TODO: don't use get()
    if (expr.get() == NULL)
    {
      expr = DiscreteTraits<discrete_object_tag>::indexed_expr_t::construct(parent, indexExpr);
    }

    return holder_t(expr);
  }
};

}

}

#endif
