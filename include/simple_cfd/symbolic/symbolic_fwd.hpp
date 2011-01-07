#ifndef SIMPLE_CFD_SYMBOLIC_SYMBOLIC_FWD_HPP
#define SIMPLE_CFD_SYMBOLIC_SYMBOLIC_FWD_HPP

#include "ertti.hpp"

#define DECLARE_SYMBOLIC_NODE(classname, parentname) \
class classname : public ERTTI< classname, parentname >, public parentname

namespace cfd
{

namespace symbolic
{

class Expr;
class Basic;
template<typename T, typename P> class ERRRI;
class Symbol;
class Number;

}

}

#endif
