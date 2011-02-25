#ifndef SIMPLE_CFD_SYMBOLIC_FLOAT_HPP
#define SIMPLE_CFD_SYMBOLIC_FLOAT_HPP

#include <string>
#include <cstddef>
#include <ostream>
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
              boost::equality_comparable<Float>
{
private:
  double value;

public:
  static Float fromFraction(long numerator, long denominator);

  Float(const double _value);
  virtual std::size_t nops() const;
  virtual void write(std::ostream& o) const;
  virtual Expr derivative(const Symbol& s) const;
  virtual bool has(const Expr& e) const;
  Expr subs(const Expr::subst_map& map) const;
  Expr integrate(const Symbol& s) const;
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
};

Float pow(const Float& f, int exponent);

}

}

#endif
