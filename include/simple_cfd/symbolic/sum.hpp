#ifndef SIMPLE_CFD_SUM_HPP
#define SIMPLE_CFD_SUM_HPP

#include <cstddef>
#include <map>
#include <ostream>
#include <utility>
#include <vector>
#include <utility>
#include <boost/operators.hpp>
#include "pair_seq.hpp"
#include "expr.hpp"
#include "rational.hpp"

namespace cfd
{

namespace symbolic
{

class Sum : public PairSeq<Sum, Rational>,
                   boost::addable<Sum>
{
protected:
  friend class PairSeq<Sum, Rational>;

  Sum(const Rational& _overall, const LazyTermMap& _terms): PairSeq<Sum, Rational>(_overall, _terms)
  {
  }

  static Rational null();
  static void combineOverall(Rational& overall, const Rational& other);
  static Rational applyCoefficient(const Rational& t, const Rational& coefficient);
  Sum extractMultipliers() const;
  
public:
  static Sum constant(const Rational& r);
  static Sum sub(const Expr& a, const Expr& b);
  static Sum rational_multiple(const Expr& e, const Rational& n);
  static Sum add(const Expr& a, const Expr& b);

  Sum()
  {
  }
  
  explicit Sum(const Expr& a);
  Rational findMultiplier() const;
  Sum& operator+=(const Expr& e);
  Sum& operator+=(const Sum& s);
  Sum& operator/=(const Rational& r);
  void write(std::ostream& o) const;
  Expr derivative(const Symbol& s, Expr::optional_derivative_cache cache) const;
  Expr integrate(const Symbol& s, unsigned flags) const;
  Sum expandedProduct(const Sum& s) const;
  Float eval(const Expr::subst_map& map) const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  Expr extractMultiplier(Rational& coeff) const;
  Expr integrate(const Expr::region_t& region, const unsigned flags) const;
};

}

}

#endif
