#ifndef SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_VALUE_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_VALUE_HPP

#include "scalar_expr.hpp"
#include <cassert>
#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace detail
{

class TemporalIndexValue
{
private:
  bool terminationSet;
  ScalarExpr::expr_ptr termination;

public:
  typedef boost::shared_ptr<TemporalIndexValue> index_ptr;

  TemporalIndexValue() : terminationSet(false)
  {
  }

  void setTermination(const ScalarExpr::expr_ptr c)
  {
    assert(!terminationSet);
    terminationSet = true;
    termination = c;
  }

  ScalarExpr& getTermination() const
  {
    return *termination;
  }
};

}

}

#endif
