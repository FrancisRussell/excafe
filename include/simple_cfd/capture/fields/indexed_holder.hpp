#ifndef SIMPLE_CFD_CAPTURE_FIELDS_INDEXED_HOLDER_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_INDEXED_HOLDER_HPP

#include "discrete_traits.hpp"
#include "temporal_index.hpp"
#include "indexed_value_helper.hpp"
#include "indexable_value.hpp"

namespace cfd
{

namespace detail
{

template<typename discrete_object_tag>
class IndexedHolder
{
private:
  typedef IndexedValueHelper<discrete_object_tag> helper_t;

  typename IndexableValue<discrete_object_tag>::value_ptr indexableValue;
  const TemporalIndex indexVariable;

public:
  IndexedHolder(const TemporalIndex& i) : 
    indexableValue(new IndexableValue<discrete_object_tag>()), indexVariable(i)
  {
  }

  helper_t operator[](unsigned int absoluteIndex)
  {
    return helper_t(indexableValue, TemporalIndexExpr::absolute(indexVariable.getIndex(), absoluteIndex));
  }

  helper_t operator[](const TemporalIndexExpr& indexExpr)
  {
    assert(indexExpr.getIndex() == indexVariable.getIndex());
    return helper_t(indexableValue, indexExpr);
  }

  helper_t operator[](const TemporalIndexOffset& offset)
  {
    const TemporalIndexExpr indexExpr(indexVariable.getIndex(), offset);
    return helper_t(indexableValue, indexExpr);
  }
};

}

typedef detail::IndexedHolder<detail::discrete_scalar_tag> IndexedScalar;
typedef detail::IndexedHolder<detail::discrete_field_tag> IndexedField;
typedef detail::IndexedHolder<detail::discrete_operator_tag> IndexedOperator;

}

#endif
