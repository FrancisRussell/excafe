#ifndef SIMPLE_CFD_SYMBOLIC_RATIONAL_HPP
#define SIMPLE_CFD_SYMBOLIC_RATIONAL_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"
#include "float.hpp"
#include <boost/operators.hpp>

namespace cfd
{

namespace symbolic
{

class Rational : public AbstractBasic<Rational>, 
                 boost::arithmetic<Rational>,
                 boost::equality_comparable<Rational>
{
private:
  long numerator;
  long denominator;

  void normalise();
  static unsigned long gcd(unsigned long a, unsigned long b);

public:
  Rational(long value);
  Rational(long numerator, long denominator);
  virtual std::size_t nops() const;
  virtual void write(std::ostream& o) const;
  virtual Expr derivative(const Symbol& s) const;
  virtual bool has(const Expr& e) const;
  Expr subs(const Expr::subst_map& map) const;
  Expr integrate(const Symbol& s) const;
  Float eval(const Expr::subst_map& map) const;
  Float toFloat() const;
  bool operator==(const Rational& n) const;
  bool operator<(const Rational& n) const;
  Rational operator-() const;
  Rational& operator+=(const Rational& r);
  Rational& operator-=(const Rational& r);
  Rational& operator/=(const Rational& r);
  Rational& operator*=(const Rational& r);
  std::size_t untypedHash() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
};

}

}

#endif
