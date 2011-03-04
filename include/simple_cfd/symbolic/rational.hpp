#ifndef SIMPLE_CFD_SYMBOLIC_RATIONAL_HPP
#define SIMPLE_CFD_SYMBOLIC_RATIONAL_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <utility>
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
public:
  typedef long value_type;

private:
  value_type numerator;
  value_type denominator;

  void normalise();
  static value_type gcd(value_type a, value_type b);
  static value_type lcm(value_type a, value_type b);

public:
  static Rational gcd(const Rational& a, const Rational& b);

  Rational();
  Rational(value_type value);
  Rational(value_type numerator, value_type denominator);
  virtual std::size_t nops() const;
  virtual void write(std::ostream& o) const;
  virtual Expr derivative(const Symbol& s) const;
  virtual bool has(const Expr& e) const;
  Expr subs(const Expr::subst_map& map) const;
  Expr integrate(const Symbol& s) const;
  Float eval(const Expr::subst_map& map) const;
  Float toFloat() const;
  value_type getNumerator() const;
  value_type getDenominator() const;
  bool operator==(const Rational& n) const;
  bool operator<(const Rational& n) const;
  Rational reciprocal() const;
  Rational operator-() const;
  Rational& operator+=(const Rational& r);
  Rational& operator-=(const Rational& r);
  Rational& operator/=(const Rational& r);
  Rational& operator*=(const Rational& r);
  Rational& operator++();
  Rational& operator--();
  std::size_t untypedHash() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
};

Rational abs(const Rational& r);
Rational pow(const Rational& r, int exponent);

}

}

#endif
