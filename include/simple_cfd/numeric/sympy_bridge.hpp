#ifndef SIMPLE_CFD_NUMERIC_SYMPY_BRIDGE_HPP
#define SIMPLE_CFD_NUMERIC_SYMPY_BRIDGE_HPP

#include <boost/python/object.hpp>

namespace cfd
{
namespace detail
{
namespace sympy_bridge
{

void init(boost::python::object& global);

enum OperatorType
{
  ADD,
  MUL,
  EXP,
  FLOAT,
  INTEGER,
  SYM,
  ABS
};

}
}
}

#endif
