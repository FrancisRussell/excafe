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

public:
  static TemporalIndexExpr absolute(const TemporalIndexValue::index_ptr& _index, const signed _offset)
  {
    return TemporalIndexExpr(_index, TemporalIndexOffset(TemporalIndexOffset::absolute_tag(), _offset));
  }

  static TemporalIndexExpr relative(const TemporalIndexValue::index_ptr& _index, const signed _offset)
  {
    return TemporalIndexExpr(_index, TemporalIndexOffset(TemporalIndexOffset::relative_tag(), _offset));
  }

  static TemporalIndexExpr final(const TemporalIndexValue::index_ptr& _index, const signed _offset)
  {
    return TemporalIndexExpr(_index, TemporalIndexOffset(TemporalIndexOffset::final_tag(), _offset));
  }

  TemporalIndexExpr(const TemporalIndexValue::index_ptr& _index, const TemporalIndexOffset& _offset) :
    index(_index), offset(_offset)
  {
  }

  TemporalIndexValue::index_ptr getIndex() const
  {
    return index;
  }

  TemporalIndexOffset getOffset() const
  {
    return offset;
  }

  bool operator==(const TemporalIndexExpr& t) const
  {
    return index == t.index && offset == t.offset;
  }

  bool operator<(const TemporalIndexExpr& t) const
  {
    if (index < t.index) return true;
    if (index == t.index && offset < t.offset) return true;
    return false;
  }
};

}

}

#endif
