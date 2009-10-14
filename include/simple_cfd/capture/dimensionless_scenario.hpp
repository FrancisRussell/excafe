#ifndef SIMPLE_CFD_CAPTURE_DIMENSIONLESS_SCENARIO_HPP
#define SIMPLE_CFD_CAPTURE_DIMENSIONLESS_SCENARIO_HPP

#include "capture_fwd.hpp"

namespace cfd
{

namespace detail
{

class DimensionlessScenario
{
public:
  virtual void execute(SolveOperation& o) = 0;
  virtual ~DimensionlessScenario() {}
};

}

}

#endif
