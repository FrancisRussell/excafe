#ifndef SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_EXPR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_EXPR_HPP

#include "temporal_index_value.hpp"
#include "temporal_index_offset.hpp"

namespace cfd
{

namespace detail
{

class TemporalIndexExpr
{
private:
  const TemporalIndexValue::index_ptr index;
  const TemporalIndexOffset offset;

  TemporalIndexExpr(const TemporalIndexValue::index_ptr& _index, const TemporalIndexOffset& _offset) :
    index(_index), offset(_offset)
  {
  }

public:
  static TemporalIndexExpr absolute(const TemporalIndexValue::index_ptr& _index, const unsigned offset)
  {
    return TemporalIndexExpr(_index, TemporalIndexOffset(TemporalIndexOffset::absolute(offset)));
  }

  static TemporalIndexExpr relative(const TemporalIndexValue::index_ptr& _index, const unsigned offset)
  {
    return TemporalIndexExpr(_index, TemporalIndexOffset(TemporalIndexOffset::relative(offset)));
  }

  static TemporalIndexExpr final(const TemporalIndexValue::index_ptr& _index)
  {
    return TemporalIndexExpr(_index, TemporalIndexOffset(TemporalIndexOffset::final()));
  }
};

}

}

#endif
