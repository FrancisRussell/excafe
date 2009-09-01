#ifndef SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_HPP

#include "temporal_index_value.hpp"
#include "temporal_index_offset.hpp"

namespace cfd
{

class TemporalIndex
{
private:
  detail::TemporalIndexValue::index_ptr value;

public:
  TemporalIndex() : value(new detail::TemporalIndexValue())
  {
  }

  detail::TemporalIndexOffset operator-(const unsigned offset) const
  {
    return detail::TemporalIndexOffset(value, offset);
  }
};

}

#endif
