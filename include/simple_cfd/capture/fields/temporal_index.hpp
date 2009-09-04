#ifndef SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_HPP

#include "temporal_index_value.hpp"
#include "temporal_index_expr.hpp"

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

  operator detail::TemporalIndexExpr() const
  {
    return detail::TemporalIndexExpr::absolute(value, 0);
  }

  detail::TemporalIndexValue::index_ptr getIndex() const
  {
    return value;
  }

  void setTermination(const Scalar& s)
  {
    value->setTermination(s.getExpr());
  }
};

}

#endif
