#ifndef EXCAFE_SYMBOLIC_FLOAT_HPP
#define EXCAFE_SYMBOLIC_FLOAT_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <set>
#include <utility>
#include <boost/operators.hpp>
#include <excafe/mp/rational.hpp>
#include <excafe/mp/float.hpp>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"

namespace excafe
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
  Float();
  Float(const double _value);
  Float(const mp::Rational& value);
  Float(const mp::Float& value);
  std::size_t nops() const;
  void write(std::ostream& o) const;
  Expr derivative(const Symbol& s, Expr::optional_derivative_cache cache) const;
  bool depends(const std::set<Symbol>& symbols) const;
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
  Float abs() const;
  std::size_t untypedHash() const;
  double toDouble() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  Expr extractMultiplier(Rational& coeff) const;
  Expr extractPolynomials(ExtractedExpressions& extracted) const;
};

Float pow(const Float& f, int exponent);

}

}

#endif
