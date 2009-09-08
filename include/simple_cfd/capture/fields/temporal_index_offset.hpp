#ifndef SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_OFFSET_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_OFFSET_HPP

#include <simple_cfd/util.hpp>
#include <boost/variant.hpp>

namespace cfd
{

namespace detail
{

class TemporalIndexOffset
{
public:
  struct absolute_tag : public cfd::util::tag<absolute_tag> {};
  struct relative_tag : public cfd::util::tag<relative_tag> {};
  struct final_tag : public cfd::util::tag<final_tag> {};

  typedef boost::variant<absolute_tag, relative_tag, final_tag> offset_t;
  offset_t offsetType;
  signed offset;

  TemporalIndexOffset(const offset_t& _offsetType, const signed _offset) : offsetType(_offsetType),
    offset(_offset)
  {
  }

  offset_t getType() const
  {
    return offsetType;
  }

  signed getValue() const
  {
    return offset;
  }

  bool operator==(const TemporalIndexOffset& t) const
  {
    return offsetType == t.offsetType && offset == t.offset;
  }

  bool operator<(const TemporalIndexOffset& t) const
  {
    if (offsetType < t.offsetType) return true;
    if (offsetType == t.offsetType && offset < t.offset) return true;
    return false;
  }
};

}

}

#endif
