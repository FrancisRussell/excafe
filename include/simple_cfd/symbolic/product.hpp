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

  static Expr integrate(const Expr& a, const Expr& b, const Symbol& s);

  Product(const Rational& _overall, const TermMap& _terms): 
    PairSeq<Product, int>(_overall, _terms)
  {
  }

  static Rational null();
  static void combineOverall(Rational& overall, const Rational& other);
  static Rational applyCoefficient(const Rational& t, const int coefficient);
  static void extractMultipliers(Rational& overall, TermMap& map);

public:
  static Product pow(const Expr& base, const int exponent)
  {
    TermMap terms;
    terms[base] += exponent;
    return Product(null(), terms);
  }

  static Product div(const Expr& a, const Expr& b)
  {
    TermMap terms;
    ++terms[a];
    --terms[b];
    return Product(null(), terms);
  }

  static Product mul(const Expr& a, const Expr& b)
  {
    TermMap terms;
    ++terms[a];
    ++terms[b];
    return Product(null(), terms);
  }

  Product& operator*=(const Product& p);
  void write(std::ostream& o) const;
  Expr derivative(const Symbol& s) const;
  Expr integrate(const Symbol& s) const;
  Expr simplify() const;
  Expr expand() const;
  Float eval(const Expr::subst_map& map) const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  virtual Expr extractMultiplier(Rational& coeff) const;
};

}

}

#endif
