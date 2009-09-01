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
  class absolute
  { 
  private:
    unsigned int offsetValue;

  public: 
    absolute(const unsigned int _offsetValue) : offsetValue(_offsetValue)
    {
    }

    unsigned int getOffsetValue() const
    {
      return offsetValue;
    }
  };

  class relative
  { 
  private:
    unsigned int offsetValue;

  public: 
    relative(const unsigned int _offsetValue) : offsetValue(_offsetValue)
    {
    }

    unsigned int getOffsetValue() const
    {
      return offsetValue;
    }
  };

  class final {};

  typedef boost::variant<absolute, relative, final> offset_t;
  offset_t offset_detail;

  TemporalIndexOffset(const offset_t& _offset_detail) : offset_detail(_offset_detail)
  {
  }

  offset_t getOffset() const
  {
    return offset_detail;
  }
};

}

}

#endif
