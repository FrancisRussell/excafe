#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_HPP

#include "function_space.hpp"

namespace cfd
{

class Operator
{
private:
  FunctionSpace trialSpace;
  FunctionSpace testSpace;

public:
  Operator(const FunctionSpace& _trialSpace, const FunctionSpace& _testSpace) : trialSpace(_trialSpace),
    testSpace(_testSpace)
  {
  }

  Operator& operator+=(int)
  {
    return *this;
  }
};

}

#endif
