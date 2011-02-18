#ifndef SIMPLE_CFD_SUM_HPP
#define SIMPLE_CFD_SUM_HPP

#include <cstddef>
#include <map>
#include <ostream>
#include <utility>
#include <vector>
#include "pair_seq.hpp"
#include "expr.hpp"

namespace cfd
{

namespace symbolic
{

class Sum : public PairSeq<Sum>
{
protected:
  friend class PairSeq<Sum>;

  Sum(const Rational& _overall, const TermMap& _terms): PairSeq<Sum>(_overall, _terms)
  {
  }

  static Rational null();
  static void combineOverall(Rational& overall, const Rational& other);
  static Rational applyCoefficient(const Rational& t, const int coefficient);
  
public:
  static Sum sub(const Expr& a, const Expr& b)
  {
    TermMap terms;
    ++terms[a];
    --terms[b];
    return Sum(null(), terms);
  }

  static Sum integer_multiple(const Expr& e, const int n)
  {
    TermMap terms;
    terms[e]+=n;
    return Sum(null(), terms);
  }

  static Sum add(const Expr& a, const Expr& b)
  {
    TermMap terms;
    ++terms[a];
    ++terms[b];
    return Sum(null(), terms);
  }

  Sum()
  {
  }
  
  explicit Sum(const Expr& a)
  {
    ++terms[a];
  }

  Sum operator+(const Expr& e) const;
  void write(std::ostream& o) const;
  Expr derivative(const Symbol& s) const;
  Expr integrate(const Symbol& s) const;
  Sum expandedProduct(const Sum& s) const;
  Expr eval() const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
};

}

}

#endif
