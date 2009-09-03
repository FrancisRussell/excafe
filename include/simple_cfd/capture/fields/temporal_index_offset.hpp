#ifndef SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_OFFSET_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_OFFSET_HPP

#include <boost/variant.hpp>

namespace cfd
{

namespace detail
{

class TemporalIndexOffset
{
public:
  struct absolute_tag {};
  struct relative_tag {};
  struct final_tag {};

  typedef boost::variant<absolute_tag, relative_tag, final_tag> offset_t;
  offset_t offsetType;
  unsigned int offset;

  TemporalIndexOffset(const offset_t& _offsetType, const unsigned _offset) : offsetType(_offsetType),
    offset(_offset)
  {
  }

  offset_t getType() const
  {
    return offsetType;
  }

  unsigned int getOffset() const
  {
    return offset;
  }
};

}

}

#endif
