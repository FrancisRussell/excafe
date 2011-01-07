#ifndef SIMPLE_CFD_SYMBOLIC_SYMBOL_HPP
#define SIMPLE_CFD_SYMBOLIC_SYMBOL_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include "basic.hpp"
#include "symbolic_fwd.hpp"
#include "expr.hpp"

namespace cfd
{

namespace symbolic
{

DECLARE_SYMBOLIC_NODE(Symbol, Basic)
{
private:
  static int nextSerial;
  std::string name;
  int serial;
  
public:
 Symbol(const std::string& _name);

 std::size_t nops() const;

 void write(std::ostream& o) const;

 Expr clone() const;

 Expr derivative(const Symbol& s) const;

 bool isNumber() const;

 std::size_t hashValue() const;
};

}

}

#endif
