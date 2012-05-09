#ifndef EXCAFE_NUMERIC_SYMPY_BRIDGE_HPP
#define EXCAFE_NUMERIC_SYMPY_BRIDGE_HPP

#include <boost/python/object.hpp>

namespace excafe
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
