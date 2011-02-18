#ifndef SIMPLE_CFD_PRODUCT_HPP
#define SIMPLE_CFD_PRODUCT_HPP

#include <cstddef>
#include <map>
#include <ostream>
#include <utility>
#include <vector>
#include "pair_seq.hpp"
#include "expr.hpp"
#include <iostream>

namespace cfd
{

namespace symbolic
{

class Product : public PairSeq<Product>
{
private:
  friend class PairSeq<Product>;

  static Expr integrate(const Expr& a, const Expr& b, const Symbol& s);

  Product(const Rational& _overall, const TermMap& _terms): 
    PairSeq<Product>(_overall, _terms)
  {
  }

  static Rational null();
  static void combineOverall(Rational& overall, const Rational& other);
  static Rational applyCoefficient(const Rational& t, const int coefficient);

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

  void write(std::ostream& o) const;
  Expr derivative(const Symbol& s) const;
  Expr integrate(const Symbol& s) const;
  Expr simplify() const;
  Expr expand() const;
  Expr eval() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
};

}

}

#endif
