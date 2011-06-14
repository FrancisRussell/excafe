#ifndef SIMPLE_CFD_SYMBOLIC_ABS_HPP
#define SIMPLE_CFD_SYMBOLIC_ABS_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <utility>
#include <set>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"

namespace cfd
{

namespace symbolic
{

class Abs : public AbstractBasic<Abs>
{
private:
  Expr expr;
  Expr sign() const;

public:
  Abs(const Expr& e);
  std::size_t nops() const;
  void write(std::ostream& o) const;
  Expr derivative(const Symbol& s, Expr::optional_derivative_cache cache) const;
  bool depends(const std::set<Symbol>& symbols) const;
  Expr getExpr() const;
  Expr subs(const Expr::subst_map& map, unsigned flags) const;
  Expr integrate(const Symbol& s, unsigned flags) const;
  Expr simplify() const;
  Float eval(const Expr::subst_map& map) const;
  bool operator==(const Abs& g) const;
  std::size_t untypedHash() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  Expr extractPolynomials(ExtractedExpressions& extracted) const;
};

}

}

#endif
