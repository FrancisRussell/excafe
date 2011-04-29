#ifndef SIMPLE_CFD_SYMBOLIC_FLOAT_HPP
#define SIMPLE_CFD_SYMBOLIC_FLOAT_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <set>
#include <utility>
#include <boost/operators.hpp>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"

namespace cfd
{

namespace symbolic
{

class Float : public AbstractBasic<Float>, 
              boost::arithmetic<Float>,
              boost::totally_ordered<Float>
{
private:
  double value;

  bool asRational(Rational& r) const;

public:
  static Float fromFraction(long numerator, long denominator);

  Float();
  Float(const double _value);
  virtual std::size_t nops() const;
  virtual void write(std::ostream& o) const;
  virtual Expr derivative(const Symbol& s) const;
  virtual bool depends(const std::set<Symbol>& symbols) const;
  Expr subs(const Expr::subst_map& map, unsigned flags) const;
  Expr integrate(const Symbol& s, unsigned flags) const;
  Expr simplify() const;
  Float eval(const Expr::subst_map& map) const;
  bool operator==(const Float& n) const;
  bool operator<(const Float& n) const;
  Float& operator+=(const Float& n);
  Float& operator-=(const Float& n);
  Float& operator*=(const Float& n);
  Float& operator/=(const Float& n);
  Float pow(int exponent) const;
  std::size_t untypedHash() const;
  double toDouble() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  virtual Expr extractMultiplier(Rational& coeff) const;
};

Float pow(const Float& f, int exponent);

}

}

#endif
