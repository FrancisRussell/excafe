#ifndef SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_OFFSET_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_OFFSET_HPP

namespace cfd
{

namespace detail
{

class TemporalIndexOffset
{
private:
  const TemporalIndexValue::index_ptr index;
  const int offset;

public:
  TemporalIndexOffset(const TemporalIndexValue::index_ptr& _index, const int _offset) :
    index(_index), offset(_offset)
  {
  }
};

}

}

#endif
