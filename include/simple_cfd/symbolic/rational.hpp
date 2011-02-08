#ifndef SIMPLE_CFD_SYMBOLIC_RATIONAL_HPP
#define SIMPLE_CFD_SYMBOLIC_RATIONAL_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"

namespace cfd
{

namespace symbolic
{

class Rational : public AbstractBasic<Rational>
{
private:
  int numerator;
  int denominator;

public:
  Rational(int value);
  Rational(int numerator, int denominator);
  virtual std::size_t nops() const;
  virtual void write(std::ostream& o) const;
  virtual Expr derivative(const Symbol& s) const;
  virtual bool has(const Expr& e) const;
  Expr subs(const Expr::subst_map& map) const;
  Expr integrate(const Symbol& s) const;
  Expr eval() const;
  bool operator==(const Rational& n) const;
  bool operator<(const Rational& n) const;
  std::size_t untypedHash() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
};

}

}

#endif
