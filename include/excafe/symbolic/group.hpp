#ifndef EXCAFE_SYMBOLIC_GROUP_HPP
#define EXCAFE_SYMBOLIC_GROUP_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <utility>
#include <set>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"

namespace excafe
{

namespace symbolic
{

class Group : public AbstractBasic<Group>
{
private:
  Expr expr;

public:
  Group(const Expr& e);
  std::size_t nops() const;
  void write(std::ostream& o) const;
  Expr derivative(const Symbol& s, Expr::optional_derivative_cache cache) const;
  bool depends(const std::set<Symbol>& symbols) const;
  Expr getExpr() const;
  Expr subs(const Expr::subst_map& map, unsigned flags) const;
  Expr integrate(const Symbol& s, unsigned flags) const;
  Expr simplify() const;
  Float eval(const Expr::subst_map& map) const;
  bool operator==(const Group& g) const;
  std::size_t untypedHash() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  bool isPolynomial() const;
  Expr extractPolynomials(ExtractedExpressions& extracted) const;
};

}

}

#endif
