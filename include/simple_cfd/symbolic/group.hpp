#ifndef SIMPLE_CFD_SYMBOLIC_GROUP_HPP
#define SIMPLE_CFD_SYMBOLIC_GROUP_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <utility>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"

namespace cfd
{

namespace symbolic
{

class Group : public AbstractBasic<Group>
{
private:
  Expr expr;

public:
  Group(const Expr& e);
  virtual std::size_t nops() const;
  virtual void write(std::ostream& o) const;
  virtual Expr derivative(const Symbol& s) const;
  virtual bool has(const Expr& e) const;
  Expr getExpr() const;
  Expr subs(const Expr::subst_map& map) const;
  Expr integrate(const Symbol& s) const;
  Expr simplify() const;
  Float eval(const Expr::subst_map& map) const;
  bool operator==(const Group& g) const;
  std::size_t untypedHash() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
};

}

}

#endif
