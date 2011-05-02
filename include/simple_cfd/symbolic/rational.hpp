#ifndef SIMPLE_CFD_SYMBOLIC_RATIONAL_HPP
#define SIMPLE_CFD_SYMBOLIC_RATIONAL_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <utility>
#include <set>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"
#include "float.hpp"
#include <boost/operators.hpp>
#include <cln/rational.h>

namespace cfd
{

namespace symbolic
{

class Rational : public AbstractBasic<Rational>, 
                 boost::arithmetic<Rational>,
                 boost::totally_ordered<Rational>
{
public:
  static const Expr zero();
  static const Expr one();

private:
  cln::cl_RA value;
  void normalise();

  Rational(const cln::cl_RA& value);
  cln::cl_I getNumerator() const;
  cln::cl_I getDenominator() const;

public:
  static Rational gcd(const Rational& a, const Rational& b);

  Rational();
  Rational(long value);
  Rational(long numerator, long denominator);
  virtual std::size_t nops() const;
  virtual void write(std::ostream& o) const;
  virtual Expr derivative(const Symbol& s) const;
  virtual bool depends(const std::set<Symbol>& symbols) const;
  Expr subs(const Expr::subst_map& map, unsigned flags) const;
  Expr integrate(const Symbol& s, unsigned flags) const;
  Float eval(const Expr::subst_map& map) const;
  Float toFloat() const;
  bool operator==(long n) const;
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
  virtual Expr extractMultiplier(Rational& coeff) const;
  Rational abs() const;
  Rational pow(int exponent) const;
};

Rational abs(const Rational& r);
Rational pow(const Rational& r, int exponent);

}

}

#endif
