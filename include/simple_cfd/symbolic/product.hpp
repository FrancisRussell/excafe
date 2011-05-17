#ifndef SIMPLE_CFD_PRODUCT_HPP
#define SIMPLE_CFD_PRODUCT_HPP

#include <cstddef>
#include <map>
#include <ostream>
#include <utility>
#include <vector>
#include <boost/operators.hpp>
#include "pair_seq.hpp"
#include "expr.hpp"
#include <iostream>

namespace cfd
{

namespace symbolic
{

class Product : public PairSeq<Product, int>,
                       boost::multipliable<Product>
{
private:
  friend class PairSeq<Product, int>;

  static Expr integrateComplex(const LazyTermMap& terms, const Symbol& s, unsigned flags);
  static Expr integrate(const Product& a, const Product& b, const Symbol& s, unsigned flags);

  Product(const Rational& _overall, const LazyTermMap& _terms): 
    PairSeq<Product, int>(_overall, _terms)
  {
  }

  static Rational null();
  static void combineOverall(Rational& overall, const Rational& other);
  static Rational applyCoefficient(const Rational& t, const int coefficient);
  Product extractMultipliers() const;

public:
  static Product constant(const Rational& r);
  static Product pow(const Expr& base, const int exponent);
  static Product div(const Expr& a, const Expr& b);
  static Product mul(const Expr& a, const Expr& b);

  Product& operator*=(const Product& p);
  void write(std::ostream& o) const;
  Expr derivative(const Symbol& s, Expr::optional_derivative_cache cache) const;
  Expr integrate(const Symbol& s, unsigned flags) const;
  Expr expand() const;
  Float eval(const Expr::subst_map& map) const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  Expr extractMultiplier(Rational& coeff) const;
  Expr integrate(const Expr::region_t& region, const unsigned flags) const;
};

}

}

#endif
