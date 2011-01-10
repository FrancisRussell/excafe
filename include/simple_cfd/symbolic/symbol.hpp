#ifndef SIMPLE_CFD_SYMBOLIC_SYMBOL_HPP
#define SIMPLE_CFD_SYMBOLIC_SYMBOL_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"
#include "basic.hpp"
#include "expr.hpp"

namespace cfd
{

namespace symbolic
{

class Symbol : public AbstractBasic<Symbol>
{
private:
  static int nextSerial;
  std::string name;
  int serial;
  
public:
 Symbol(const std::string& _name);
 std::size_t nops() const;
 void write(std::ostream& o) const;
 Expr derivative(const Symbol& s) const;
 bool isNumber() const;
 bool operator==(const Symbol& s) const;
 bool operator<(const Symbol& s) const;
 Expr integrate(const Symbol& s) const;
 Expr subs(const Expr::subst_map& map) const;
 std::size_t untypedHash() const;
};

}

}

#endif
