#ifndef SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_VALUE_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_VALUE_HPP

#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace detail
{

class TemporalIndexValue
{
public:
  typedef boost::shared_ptr<TemporalIndexValue> index_ptr;

  TemporalIndexValue()
  {
  }
};

}

}

#endif
