#ifndef SIMPLE_CFD_SYMBOLIC_SYMBOL_HPP
#define SIMPLE_CFD_SYMBOLIC_SYMBOL_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <utility>
#include <set>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"
#include "basic.hpp"
#include "expr.hpp"
#include "rational.hpp"
#include <boost/operators.hpp>

namespace cfd
{

namespace symbolic
{

class Symbol : public AbstractBasic<Symbol>,
               boost::totally_ordered<Symbol>
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
  Float eval(const Expr::subst_map& map) const;
  bool operator==(const Symbol& s) const;
  bool operator<(const Symbol& s) const;
  bool depends(const std::set<Symbol>& symbols) const;
  Expr integrate(const Symbol& s, unsigned flags) const;
  Expr subs(const Expr::subst_map& map, unsigned flags) const;
  std::size_t untypedHash() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
};

}

}

#endif
